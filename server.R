server <- function(input, output, session) {

  # ── Landing page "Start Analysis" button ─────────────────────────────────
  observeEvent(input$go_to_home, {
    updateNavbarPage(session, "main_nav", selected = "Analysis")
    shinyjs::runjs("window.scrollTo(0, 0);")
  })

  # ── About page "How It Works" link ───────────────────────────────────────
  observeEvent(input$about_to_hiw, {
    updateNavbarPage(session, "main_nav", selected = "How It Works")
    shinyjs::runjs("window.scrollTo(0, 0);")
  })

  # ════════════════════════════════════════════════════════════════════════════
  #  "How It Works" toy examples — plain arithmetic only, no Stan/MalReBay
  #  calls, but mirroring the real model's math (src/stan/malrebay_model.stan,
  #  dist_log_prob_logic()): method 1 (microsatellite / length-polymorphic)
  #  uses a distance-decay kernel; method 3 (AmpSeq / microhaplotype) uses an
  #  exact-match-only genotyping-error model. Both share the same downstream
  #  structure: LR = P(observation | same lineage) / freq(recurrence allele),
  #  combined across Day-0 x recurrence allele pairs via a log-mean-exp
  #  exactly matching log_sum_exp(lpr) - log(MOI0*MOIf) in the .stan model.
  # ════════════════════════════════════════════════════════════════════════════
  hiw_K <- 10  # number of possible alleles/haplotypes for the toy marker

  # Shared: a geometric-shaped population frequency distribution over K
  # alleles, controlled by "evenness" (higher = flatter / more diverse).
  hiw_freq_from_evenness <- function(evenness) {
    skew <- 1 - evenness
    raw  <- (1 - skew)^(seq_len(hiw_K) - 1)
    raw / sum(raw)
  }

  # Shared: combine per-pair likelihood ratios into one marker-level value.
  hiw_log_mean_exp_lr <- function(lrs) {
    log_lrs <- log(lrs)
    n <- length(log_lrs)
    m <- max(log_lrs)
    exp(m + log(sum(exp(log_lrs - m))) - log(n))
  }

  # Shared: an inline-styled badge (not a Bootstrap badge-* class, which
  # doesn't render visibly against this app's theme) in the same red/blue
  # color language used for Recrudescence/New Infection elsewhere in the app.
  hiw_badge <- function(text, is_red) {
    span(style = paste0(
           "display:inline-block; padding:0.2em 0.6em; border-radius:0.25rem;",
           "font-weight:700; font-size:0.85em; color:#fff; background-color:",
           if (is_red) "#D6604D" else "#67A9CF", ";"),
         text)
  }

  # Shared: the 3-marker match-counting comparison table (strict 3/3 vs
  # loose 2/3), given marker 1's derived status and the two toggle markers.
  hiw_match_counting_table <- function(m1, m2, m3) {
    match_count <- sum(m1, m2, m3)

    status_badge <- function(is_match) hiw_badge(if (is_match) "Match" else "Mismatch", is_match)
    result_badge <- function(is_r) hiw_badge(if (is_r) "Recrudescence" else "New infection", is_r)

    tagList(
      div(class = "mb-3",
        tags$strong("Marker 1: "), status_badge(m1), "  ",
        tags$strong("Marker 2: "), status_badge(m2), "  ",
        tags$strong("Marker 3: "), status_badge(m3)
      ),
      tags$table(class = "table table-sm table-bordered",
        tags$thead(tags$tr(
          tags$th("Rule"), tags$th("Threshold"), tags$th("Markers matched"), tags$th("Classification")
        )),
        tags$tbody(
          tags$tr(
            tags$td("Strict"), tags$td("3/3 must match"),
            tags$td(paste0(match_count, "/3")),
            tags$td(result_badge(match_count == 3))
          ),
          tags$tr(
            tags$td("Loose"), tags$td("2/3 must match"),
            tags$td(paste0(match_count, "/3")),
            tags$td(result_badge(match_count >= 2))
          )
        )
      )
    )
  }

  # Shared: step 3 of the result panel - the marker's combined strength of
  # evidence, as a number and a sentence, with the per-pair calculation behind a disclosure.
  # `pairs_display` holds the table columns to show; `pair_lrs` the matching
  # likelihood ratios.
  hiw_evidence_ui <- function(pairs_display, pair_lrs, marker_lr, allele_word) {
    # Headline: the factor in the favoured direction, so the number and the
    # sentence agree (an LR of 0.40 is shown as 2.5× for a new infection)
    favours_recrud <- marker_lr >= 1
    factor   <- if (favours_recrud) marker_lr else 1 / marker_lr
    sentence <- sprintf("The data are %.1f× more likely under %s.", factor,
                        if (favours_recrud) "recrudescence" else "a new infection")
    combined_note <- sprintf(
      "Averaging the pairs gives LR = %s for the marker. %s",
      signif(marker_lr, 3),
      if (favours_recrud) "Above 1 favours recrudescence; this is the number shown above."
      else "Below 1 favours a new infection; the number shown above is 1 ÷ LR.")

    header <- tags$tr(lapply(c(names(pairs_display), "LR"), tags$th))
    rows <- lapply(seq_len(nrow(pairs_display)), function(i) {
      tags$tr(
        lapply(pairs_display[i, ], function(v) tags$td(as.character(v))),
        tags$td(sprintf("%.2f×", pair_lrs[i]))
      )
    })

    tagList(
      div(class = "hiw-subhead", span(class = "hiw-step-badge", "3"), " Strength of the evidence"),
      div(class = "hiw-evidence-number", sprintf("%.1f×", factor)),
      p(class = "hiw-evidence-text", sentence),
      tags$details(class = "hiw-details",
        tags$summary("Show the calculation"),
        div(class = "table-responsive mt-2",
          tags$table(class = "table table-sm table-bordered mb-0",
            tags$thead(header), tags$tbody(rows))),
        p(class = "hiw-help", HTML(paste0(
          "<strong>L(R)</strong>: likelihood of the comparison if it is a ",
          "recrudescence (the chart in step 2). <strong>L(NI)</strong>: likelihood ",
          "if it is a new infection. The recurrence ", sub("s$", "", allele_word),
          " is then a random draw from the local parasites, so this is its local ",
          "frequency (step 1). <strong>LR</strong> (likelihood ratio) = L(R) ÷ L(NI). ",
          combined_note)))
      )
    )
  }

  # Shared: step 4 of the result panel - starting assumption x strength of
  # evidence -> probability of recrudescence (0-1)
  hiw_post_prob <- function(prior, marker_lr) {
    prior_odds <- prior / (1 - prior)
    post_odds  <- prior_odds * marker_lr
    post_odds / (1 + post_odds)
  }

  # The probability as a gauge, with the classification threshold marked
  # (probabilities shown as 0-1 values; the CSS positions still need percent)
  hiw_prob_ui <- function(post_prob, threshold) {
    pct       <- round(post_prob * 100, 1)
    gradient  <- grDevices::colorRampPalette(c("#67A9CF", "#F7F7F7", "#D6604D"))(101)
    bar_color <- gradient[max(1, min(101, round(pct) + 1))]

    tagList(
      div(class = "hiw-subhead", span(class = "hiw-step-badge", "4"), " Probability of recrudescence"),
      div(class = "hiw-gauge-track",
        div(class = "hiw-gauge-fill",
            style = paste0("width:", pct, "%; background:", bar_color, ";")),
        div(class = "hiw-gauge-threshold", style = paste0("left:", threshold * 100, "%;"))
      ),
      div(class = "hiw-gauge-label", sprintf("%.2f", post_prob))
    )
  }

  # The classification at the chosen threshold
  hiw_verdict_ui <- function(post_prob, threshold) {
    is_recrud <- post_prob >= threshold
    tagList(
      div(class = "hiw-verdict",
        "Classification: ", hiw_badge(if (is_recrud) "Recrudescence" else "New infection", is_recrud)),
      p(class = "hiw-help", style = "text-align:center;",
        sprintf("%.2f or more \u2192 recrudescence.", threshold))
    )
  }

  # Shared: a probability/frequency for display, 2 significant figures
  hiw_fmt_prob <- function(p) ifelse(p < 0.001, "<0.001", as.character(signif(p, 2)))

  # Shared: bar-chart styling for the "How It Works" plots
  hiw_theme <- function() {
    theme_minimal(base_size = 13) +
      theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
            legend.position = "top")
  }

  # ── Length Polymorphic (microsatellite): distance-decay kernel ─────────────

  hiw_lp_freq_vec <- reactive({
    req(input$hiw_lp_evenness)
    hiw_freq_from_evenness(input$hiw_lp_evenness)
  })

  # Distance-decay kernel: P(observed size difference = m | same lineage),
  # m = 0 .. hiw_K-1, identical in form to the .stan model's dvect.
  hiw_lp_kernel_vec <- reactive({
    req(input$hiw_lp_dparam)
    d   <- input$hiw_lp_dparam
    m   <- seq(0, hiw_K - 1)
    raw <- (1 - d)^m * d
    raw / sum(raw)
  })

  hiw_lp_d0 <- reactive({
    req(length(input$hiw_lp_d0) > 0)
    as.integer(input$hiw_lp_d0)
  })
  hiw_lp_dx <- reactive({
    req(length(input$hiw_lp_dx) > 0)
    as.integer(input$hiw_lp_dx)
  })

  hiw_lp_pairs <- reactive({
    d0     <- hiw_lp_d0()
    dx     <- hiw_lp_dx()
    freq   <- hiw_lp_freq_vec()
    kernel <- hiw_lp_kernel_vec()

    grid <- expand.grid(D0 = d0, Recurrence = dx, KEEP.OUT.ATTRS = FALSE)
    grid$Distance <- abs(grid$D0 - grid$Recurrence)
    grid$P_dist   <- kernel[grid$Distance + 1]
    grid$Freq     <- freq[grid$Recurrence]
    grid$LR       <- grid$P_dist / grid$Freq
    grid
  })

  hiw_lp_marker_lr <- reactive(hiw_log_mean_exp_lr(hiw_lp_pairs()$LR))

  output$hiw_lp_freq_plot <- renderPlot({
    freq <- hiw_lp_freq_vec()
    df   <- data.frame(Allele = seq_len(hiw_K), Freq = freq)

    marks <- dplyr::bind_rows(
      data.frame(Allele = hiw_lp_d0(), Group = "Day 0"),
      data.frame(Allele = hiw_lp_dx(), Group = "Recurrence")
    )
    # Stack the two triangles when an allele is found at both timepoints
    marks$y <- freq[marks$Allele] + max(freq) * ifelse(marks$Group == "Day 0", 0.06, 0.15)

    ggplot(df, aes(x = Allele, y = Freq)) +
      geom_col(fill = "#BFC9D9", width = 0.65) +
      geom_point(data = marks, aes(x = Allele, y = y, color = Group),
                 size = 4, shape = 17, inherit.aes = FALSE) +
      scale_color_manual(values = c("Day 0" = "#2C3E50", "Recurrence" = "#D6604D")) +
      scale_x_continuous(breaks = seq_len(hiw_K)) +
      scale_y_continuous() +
      labs(x = "Allele (size)", y = "Share of local parasites", color = NULL) +
      hiw_theme()
  })

  # Labels marking the patient's own outcomes on a "likelihood if recrudescence"
  # bar chart: a marker plus "<pairs>: <chance>" above each bar the patient
  # has, so they stay visible even when the bar itself is only a sliver
  hiw_observed_labels <- function(x, p, pair_labels, y_max) {
    data.frame(x = x, P = p, label = paste0(pair_labels, "\n", hiw_fmt_prob(p)), y = p + y_max * 0.06)
  }

  # Probability axis whose headroom (for those labels) doesn't add ticks past 1
  hiw_prob_axis <- function() {
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.35)),
                       breaks = function(lims) Filter(function(b) b <= 1, scales::breaks_pretty()(lims)))
  }

  output$hiw_lp_kernel_plot <- renderPlot({
    kernel <- hiw_lp_kernel_vec()
    pairs  <- hiw_lp_pairs()
    df <- data.frame(Distance = seq(0, hiw_K - 1), P = kernel)
    df$Observed <- ifelse(df$Distance %in% pairs$Distance, "Found in this patient", "Other")

    # One label per observed difference, listing the allele pairs behind it
    pair_text <- tapply(paste(pairs$D0, "vs", pairs$Recurrence), pairs$Distance,
                        paste, collapse = "\n")
    obs_dist  <- as.integer(names(pair_text))
    marks     <- hiw_observed_labels(obs_dist, kernel[obs_dist + 1], pair_text, max(kernel))

    ggplot(df, aes(x = Distance, y = P)) +
      geom_col(aes(fill = Observed), width = 0.65) +
      geom_point(data = marks, aes(x = x, y = y), shape = 25, size = 3,
                 fill = "#D6604D", color = "#D6604D", inherit.aes = FALSE) +
      geom_text(data = marks, aes(x = x, y = y, label = label), vjust = -0.35,
                size = 3.4, lineheight = 0.9, color = "#2C3E50", inherit.aes = FALSE) +
      scale_fill_manual(values = c("Found in this patient" = "#D6604D", "Other" = "#BFC9D9"),
                        breaks = "Found in this patient", name = NULL) +
      scale_x_continuous(breaks = df$Distance) +
      hiw_prob_axis() +
      labs(x = "Difference in size between a Day 0 and a recurrence allele",
           y = "Likelihood if recrudescence, L(R)") +
      hiw_theme()
  })

  output$hiw_lp_evidence_ui <- renderUI({
    pairs <- hiw_lp_pairs()
    hiw_evidence_ui(
      data.frame(
        `Allele comparison` = paste(pairs$D0, "vs", pairs$Recurrence),
        `Size difference`   = pairs$Distance,
        `L(R)`              = hiw_fmt_prob(pairs$P_dist),
        `L(NI)`             = hiw_fmt_prob(pairs$Freq),
        check.names = FALSE
      ),
      pairs$LR, hiw_lp_marker_lr(), "alleles"
    )
  })

  hiw_lp_post_prob <- reactive({
    req(input$hiw_lp_prior)
    hiw_post_prob(input$hiw_lp_prior, hiw_lp_marker_lr())
  })

  output$hiw_lp_prob_ui <- renderUI({
    req(input$hiw_lp_threshold)
    hiw_prob_ui(hiw_lp_post_prob(), input$hiw_lp_threshold)
  })

  output$hiw_lp_verdict_ui <- renderUI({
    req(input$hiw_lp_threshold)
    hiw_verdict_ui(hiw_lp_post_prob(), input$hiw_lp_threshold)
  })

  hiw_lp_marker1_match <- reactive(any(hiw_lp_pairs()$Distance == 0))

  output$hiw_lp_match_counting_ui <- renderUI({
    req(!is.null(input$hiw_lp_match_2), !is.null(input$hiw_lp_match_3))
    hiw_match_counting_table(
      hiw_lp_marker1_match(),
      identical(input$hiw_lp_match_2, "match"),
      identical(input$hiw_lp_match_3, "match")
    )
  })

  # ── Microhaplotype / AmpSeq: exact-match-only genotyping-error model ───────

  hiw_hap_freq_vec <- reactive({
    req(input$hiw_hap_evenness)
    hiw_freq_from_evenness(input$hiw_hap_evenness)
  })

  hiw_hap_d0 <- reactive({
    req(length(input$hiw_hap_d0) > 0)
    match(input$hiw_hap_d0, paste0("H", seq_len(hiw_K)))
  })
  hiw_hap_dx <- reactive({
    req(length(input$hiw_hap_dx) > 0)
    match(input$hiw_hap_dx, paste0("H", seq_len(hiw_K)))
  })

  # Every Day-0-haplotype x recurrence-haplotype pair. A match is exact or
  # not - no distance - governed by a single genotyping/mutation error rate
  # qq: P(same | same lineage) = 1 - qq; P(this other one | same lineage) =
  # qq / (K-1), exactly matching method 3 in dist_log_prob_logic().
  hiw_hap_pairs <- reactive({
    req(input$hiw_hap_qq)
    d0   <- hiw_hap_d0()
    dx   <- hiw_hap_dx()
    freq <- hiw_hap_freq_vec()
    qq   <- input$hiw_hap_qq

    grid <- expand.grid(D0 = d0, Recurrence = dx, KEEP.OUT.ATTRS = FALSE)
    grid$Same  <- grid$D0 == grid$Recurrence
    grid$P_obs <- ifelse(grid$Same, 1 - qq, qq / (hiw_K - 1))
    grid$Freq  <- freq[grid$Recurrence]
    grid$LR    <- grid$P_obs / grid$Freq
    grid$D0_label         <- paste0("H", grid$D0)
    grid$Recurrence_label <- paste0("H", grid$Recurrence)
    grid
  })

  hiw_hap_marker_lr <- reactive(hiw_log_mean_exp_lr(hiw_hap_pairs()$LR))

  output$hiw_hap_freq_plot <- renderPlot({
    freq <- hiw_hap_freq_vec()
    df   <- data.frame(Haplotype = factor(paste0("H", seq_len(hiw_K)),
                                           levels = paste0("H", seq_len(hiw_K))),
                        Freq = freq)

    marks <- dplyr::bind_rows(
      data.frame(Haplotype = paste0("H", hiw_hap_d0()), Group = "Day 0"),
      data.frame(Haplotype = paste0("H", hiw_hap_dx()), Group = "Recurrence")
    )
    marks$Haplotype <- factor(marks$Haplotype, levels = levels(df$Haplotype))
    marks$y <- freq[match(marks$Haplotype, df$Haplotype)] +
      max(freq) * ifelse(marks$Group == "Day 0", 0.06, 0.15)

    ggplot(df, aes(x = Haplotype, y = Freq)) +
      geom_col(fill = "#BFC9D9", width = 0.65) +
      geom_point(data = marks, aes(x = Haplotype, y = y, color = Group),
                 size = 4, shape = 17, inherit.aes = FALSE) +
      scale_color_manual(values = c("Day 0" = "#2C3E50", "Recurrence" = "#D6604D")) +
      scale_y_continuous() +
      labs(x = "Haplotype", y = "Share of local parasites", color = NULL) +
      hiw_theme()
  })

  output$hiw_hap_error_plot <- renderPlot({
    req(input$hiw_hap_qq)
    qq    <- input$hiw_hap_qq
    pairs <- hiw_hap_pairs()
    outcomes <- c("Identical haplotype", "A specific different haplotype")
    df <- data.frame(
      Outcome = factor(outcomes, levels = outcomes),
      P       = c(1 - qq, qq / (hiw_K - 1))
    )
    observed <- c(if (any(pairs$Same)) outcomes[1], if (any(!pairs$Same)) outcomes[2])
    df$Observed <- ifelse(df$Outcome %in% observed, "Found in this patient", "Other")

    pair_text <- tapply(paste(pairs$D0_label, "vs", pairs$Recurrence_label),
                        outcomes[ifelse(pairs$Same, 1, 2)], paste, collapse = "\n")
    obs_idx   <- match(names(pair_text), outcomes)
    marks     <- hiw_observed_labels(obs_idx, df$P[obs_idx], pair_text, max(df$P))

    ggplot(df, aes(x = Outcome, y = P)) +
      geom_col(aes(fill = Observed), width = 0.5) +
      geom_point(data = marks, aes(x = x, y = y), shape = 25, size = 3,
                 fill = "#D6604D", color = "#D6604D", inherit.aes = FALSE) +
      geom_text(data = marks, aes(x = x, y = y, label = label), vjust = -0.35,
                size = 3.4, lineheight = 0.9, color = "#2C3E50", inherit.aes = FALSE) +
      scale_fill_manual(values = c("Found in this patient" = "#D6604D", "Other" = "#BFC9D9"),
                        breaks = "Found in this patient", name = NULL) +
      hiw_prob_axis() +
      labs(x = NULL, y = "Likelihood if recrudescence, L(R)") +
      hiw_theme()
  })

  output$hiw_hap_evidence_ui <- renderUI({
    pairs <- hiw_hap_pairs()
    hiw_evidence_ui(
      data.frame(
        `Haplotype comparison` = paste(pairs$D0_label, "vs", pairs$Recurrence_label),
        `Identical?`           = ifelse(pairs$Same, "Yes", "No"),
        `L(R)`                 = hiw_fmt_prob(pairs$P_obs),
        `L(NI)`                = hiw_fmt_prob(pairs$Freq),
        check.names = FALSE
      ),
      pairs$LR, hiw_hap_marker_lr(), "haplotypes"
    )
  })

  hiw_hap_post_prob <- reactive({
    req(input$hiw_hap_prior)
    hiw_post_prob(input$hiw_hap_prior, hiw_hap_marker_lr())
  })

  output$hiw_hap_prob_ui <- renderUI({
    req(input$hiw_hap_threshold)
    hiw_prob_ui(hiw_hap_post_prob(), input$hiw_hap_threshold)
  })

  output$hiw_hap_verdict_ui <- renderUI({
    req(input$hiw_hap_threshold)
    hiw_verdict_ui(hiw_hap_post_prob(), input$hiw_hap_threshold)
  })

  hiw_hap_marker1_match <- reactive(any(hiw_hap_pairs()$Same))

  output$hiw_hap_match_counting_ui <- renderUI({
    req(!is.null(input$hiw_hap_match_2), !is.null(input$hiw_hap_match_3))
    hiw_match_counting_table(
      hiw_hap_marker1_match(),
      identical(input$hiw_hap_match_2, "match"),
      identical(input$hiw_hap_match_3, "match")
    )
  })

  rv <- reactiveValues(
    imported_data = NULL,
    results       = NULL,
    import_error  = NULL,
    run_error     = NULL,
    example_geno   = NULL,
    example_marker = NULL,
    dataset_label  = NULL,
    run_info       = NULL   # elapsed time and MCMC settings of the last run
  )

  # ── Load example dataset helper ──────────────────────────────────────────
  load_example <- function(geno_file, marker_file = NULL, label,
                            n_chains = DEFAULT_MCMC$n_chains,
                            iter     = DEFAULT_MCMC$iter) {
    rv$import_error  <- NULL
    rv$imported_data <- NULL
    rv$results       <- NULL
    rv$example_geno   <- geno_file
    rv$example_marker <- marker_file
    rv$dataset_label  <- label

    updateNumericInput(session, "n_chains",     value = n_chains)
    updateNumericInput(session, "iter",         value = iter)
    updateNumericInput(session, "burn_in_frac", value = DEFAULT_MCMC$burn_in_frac)
    updateNumericInput(session, "adapt_delta",  value = DEFAULT_MCMC$adapt_delta)
    updateNumericInput(session, "random_seed",  value = DEFAULT_MCMC$random_seed)

    args <- list(filepath = geno_file, verbose = FALSE)
    if (!is.null(marker_file)) args$marker_filepath <- marker_file

    result <- tryCatch(do.call(import_data, args), error = function(e) e)

    if (inherits(result, "error")) {
      rv$import_error <- result$message
      showNotification(paste("Import error:", result$message),
                       type = "error", duration = 15)
    } else {
      rv$imported_data <- result
      showNotification("Example dataset loaded.", type = "message", duration = 4)
      switch_to_data()
    }
  }

  observeEvent(input$load_example_lp_short, {
    load_example(
      geno_file   = build_lp_example(LP_EXAMPLE_SHORT_SITES, "example_length_polymorphic_1site.xlsx"),
      marker_file = file.path(EXTDATA_PATH, "makers_details.xlsx"),
      label       = "Example: Length Polymorphic - 1 site",
      iter        = 1000
    )
  })

  observeEvent(input$load_example_lp, {
    load_example(
      geno_file   = build_lp_example(LP_EXAMPLE_LONG_SITES, "example_length_polymorphic_2sites.xlsx"),
      marker_file = file.path(EXTDATA_PATH, "makers_details.xlsx"),
      label       = "Example: Length Polymorphic - 2 sites",
      iter        = 1000
    )
  })

  observeEvent(input$load_example_ampseq, {
    load_example(
      geno_file = file.path(EXTDATA_PATH, "Dataset_amplicon_sequencing.xlsx"),
      label     = "Example: Amplicon Sequencing - 1 site"
    )
  })

  # ── Sidebar collapse / expand ────────────────────────────────────────────
  observeEvent(input$sidebar_toggle, {
    shinyjs::toggleClass("sidebar-col", "sidebar-collapsed")
  })

  # ── Panel toggle ─────────────────────────────────────────────────────────
  # Three mutually-exclusive main panels, switched by the toggle bar buttons.
  show_panel <- function(panel_id, btn_id) {
    all_panels <- c("data_exploration_panel", "classification_panel", "efficacy_panel")
    all_btns   <- c("btn_data_panel", "btn_results_panel", "btn_efficacy_panel")

    for (p in all_panels) if (p == panel_id) shinyjs::show(p) else shinyjs::hide(p)
    for (b in all_btns)   if (b == btn_id)   shinyjs::addClass(b, "active-panel") else shinyjs::removeClass(b, "active-panel")
  }

  observeEvent(input$btn_data_panel,     show_panel("data_exploration_panel", "btn_data_panel"))
  observeEvent(input$btn_results_panel,  show_panel("classification_panel",   "btn_results_panel"))
  observeEvent(input$btn_efficacy_panel, show_panel("efficacy_panel",         "btn_efficacy_panel"))

  switch_to_results <- function() show_panel("classification_panel",   "btn_results_panel")
  switch_to_data     <- function() show_panel("data_exploration_panel", "btn_data_panel")

  # ── Enable / disable Run button ──────────────────────────────────────────
  observe({
    files_ready   <- !is.null(input$genotype_file) && !is.null(input$marker_file)
    example_ready <- !is.null(rv$example_geno)
    if ((files_ready || example_ready) && !run_active()) {
      shinyjs::enable("run_btn")
    } else {
      shinyjs::disable("run_btn")
    }
    updateActionButton(session, "run_btn",
                       label = if (run_active()) " Running..." else " Run MalReBay",
                       icon  = if (run_active()) icon("circle-notch", class = "fa-spin") else icon("play"))
  })

  # ── Auto-import when both files are uploaded ─────────────────────────────
  observe({
    req(input$genotype_file, input$marker_file)
    rv$example_geno   <- NULL
    rv$example_marker <- NULL
    shinyjs::removeClass("load_example_lp_short", "btn-example-selected")
    shinyjs::removeClass("load_example_lp",    "btn-example-selected")
    shinyjs::removeClass("load_example_ampseq", "btn-example-selected")

    rv$import_error  <- NULL
    rv$imported_data <- NULL

    result <- tryCatch(
      import_data(
        filepath        = input$genotype_file$datapath,
        marker_filepath = input$marker_file$datapath,
        verbose         = FALSE
      ),
      error = function(e) e
    )

    if (inherits(result, "error")) {
      rv$import_error <- result$message
      showNotification(paste("Import error:", result$message),
                       type = "error", duration = 15)
    } else {
      rv$imported_data <- result
      rv$dataset_label <- input$genotype_file$name
      showNotification("Data imported successfully.", type = "message", duration = 4)
      switch_to_data()
    }
  })

  # ── Run full MalReBay pipeline ───────────────────────────────────────────
  # The run happens in a background R process (see run_malrebay_background()
  # in global.R): the app stays responsive while it runs, for this user and
  # for everyone else on the same worker.
  run_task <- ExtendedTask$new(function(run_args, started_flag)
    run_malrebay_background(run_args, started_flag))

  # The run in progress: start time, dataset, settings, and whether it had to
  # wait for a free worker. NULL when no run is going on.
  run_state <- reactiveVal(NULL)
  run_active <- reactive(!is.null(run_state()))

  # Dataset size: recurrent patients and additional Day 0 samples
  dataset_size <- reactive({
    req(rv$imported_data)
    list(patients   = length(unique(gsub(" Day 0$| recurrence$", "", rv$imported_data$late_failures$Sample.ID))),
         additional = NROW(rv$imported_data$additional))
  })

  set_result_spinners <- function(show) {
    for (id in c("spinner_patient_results", "spinner_comparison", "spinner_convergence",
                 "spinner_site_stats", "spinner_prob_histogram")) {
      if (show) shinyjs::show(id) else shinyjs::hide(id)
    }
  }

  observeEvent(input$run_btn, {
    req(rv$imported_data)
    if (run_active()) return()

    # The input's max isn't enforced on typed values (or an emptied box):
    # clamp to 1..MAX_MCMC_CHAINS here too
    clamp_chains <- function(n) {
      if (is.null(n) || is.na(n)) n <- DEFAULT_MCMC$n_chains
      min(max(1L, as.integer(round(n))), MAX_MCMC_CHAINS)
    }
    if (!is.null(input$n_chains) && !is.na(input$n_chains) && input$n_chains > MAX_MCMC_CHAINS) {
      updateNumericInput(session, "n_chains", value = MAX_MCMC_CHAINS)
      showNotification(sprintf("At most %d chains can be used; running with %d.",
                               MAX_MCMC_CHAINS, MAX_MCMC_CHAINS), type = "warning", duration = 6)
    }
    iter <- if (is.null(input$iter) || is.na(input$iter)) DEFAULT_MCMC$iter else as.integer(input$iter)

    size    <- dataset_size()
    problem <- run_size_problem(size$patients, size$additional, iter)
    if (!is.null(problem)) {
      showModal(modalDialog(
        title = "This analysis is too large for the dashboard",
        p(problem),
        p("Split the dataset by site and analyse each part separately, or use the MalReBay R package
           for larger studies; see the About page for the package website, vignette and analysis
           notebooks."),
        easyClose = TRUE, footer = modalButton("OK")))
      return()
    }

    mcmc_config <- list(
      n_chains     = clamp_chains(input$n_chains),
      iter         = iter,
      burn_in_frac = input$burn_in_frac,
      adapt_delta  = input$adapt_delta,
      random_seed  = input$random_seed
    )

    geno_path   <- if (!is.null(rv$example_geno)) rv$example_geno   else input$genotype_file$datapath
    marker_path <- if (!is.null(rv$example_marker)) rv$example_marker else input$marker_file$datapath

    # prob_threshold only seeds the package's own Classification columns;
    # the tables below re-derive them from the slider, so moving it after
    # a run doesn't require re-running the MCMC.
    # plots = FALSE: nothing is written to disk; every chart (convergence
    # included) is rebuilt from the returned results by the outputs below.
    run_args <- list(
      filepath       = normalizePath(geno_path),
      mcmc_config    = mcmc_config,
      output_folder  = NULL,
      verbose        = FALSE,
      prob_threshold = input$classification_threshold,
      plots          = FALSE
    )
    if (!is.null(marker_path)) run_args$marker_filepath <- normalizePath(marker_path)

    rv$results   <- NULL
    rv$run_error <- NULL
    started_flag <- tempfile("malrebay_run_started_")
    run_state(list(
      submitted    = Sys.time(),
      started_flag = started_flag,
      dataset      = rv$imported_data,   # to drop the results if other data are loaded meanwhile
      mcmc_config  = mcmc_config
    ))
    set_result_spinners(TRUE)
    run_task$invoke(run_args, started_flag)

    # Grey out this user's page until the run finishes (other users of the
    # app are not affected: the run happens in a background process)
    showModal(modalDialog(
      title = NULL, footer = NULL, easyClose = FALSE,
      uiOutput("run_modal_ui")))
  })

  # Run finished (or failed)
  observeEvent(run_task$status(), {
    status <- run_task$status()
    if (!status %in% c("success", "error") || !run_active()) return()
    state <- run_state()
    run_state(NULL)
    unlink(state$started_flag)
    removeModal()
    set_result_spinners(FALSE)

    outcome <- tryCatch(run_task$result(), error = function(e) e)
    if (inherits(outcome, "error")) {
      rv$run_error <- conditionMessage(outcome)
      showNotification(paste("Error:", rv$run_error), type = "error", duration = 20)
      return()
    }
    if (!identical(state$dataset, rv$imported_data)) {
      showNotification("The run finished, but another dataset has been loaded since: its results were discarded.",
                       type = "warning", duration = 10)
      return()
    }
    rv$results  <- outcome$results
    rv$run_info <- list(seconds = outcome$seconds, n_chains = state$mcmc_config$n_chains,
                        iter = state$mcmc_config$iter)
    showNotification("Analysis complete!", type = "message", duration = 5)
    switch_to_results()
  })

  # Phase of the run in progress, re-checked every second: "queued" while all
  # background workers are busy with other users' runs, then "running";
  # with the time spent in the current phase
  run_phase <- reactive({
    req(run_active())
    invalidateLater(1000)
    state <- run_state()
    started <- file.exists(state$started_flag)
    since   <- if (started) file.mtime(state$started_flag) else state$submitted
    list(phase   = if (started) "running" else "queued",
         elapsed = as.numeric(difftime(Sys.time(), since, units = "secs")))
  })

  # Content of the grey pop-up shown to the user who pressed Run
  output$run_modal_ui <- renderUI({
    req(run_active())
    ph <- run_phase()
    div(style = "text-align:center; padding:1.5rem 0.5rem;",
      icon(if (ph$phase == "queued") "hourglass-half" else "circle-notch",
           class = if (ph$phase == "running") "fa-spin",
           style = "font-size:3rem; color:#2C3E50;"),
      if (ph$phase == "queued") tagList(
        h4(style = "margin-top:1.2rem; color:#2C3E50; font-weight:700;",
           "Waiting for a free worker…"),
        p(style = "color:#7F8C8D; margin-top:0.5rem;",
          sprintf("All %d analysis workers are busy with other users' runs.", N_RUN_WORKERS),
          br(), "Your analysis is queued and will start automatically as soon as one is free.",
          br(), sprintf("Waiting for %s.", format_duration(ph$elapsed)))
      ) else tagList(
        h4(style = "margin-top:1.2rem; color:#2C3E50; font-weight:700;",
           "Running MalReBay…"),
        p(style = "color:#7F8C8D; margin-top:0.5rem;",
          "This will take a couple of minutes while the algorithm runs.",
          br(), "You can go grab a coffee ;-)"),
        p(style = "color:#95A5A6; font-size:0.9rem;", sprintf("Elapsed: %s", format_duration(ph$elapsed)))
      )
    )
  })

  # ── Sidebar: run status (under the Run button) ───────────────────────────
  output$run_status_ui <- renderUI({
    if (run_active()) {
      ph <- run_phase()
      return(div(class = "status-strip alert alert-info",
        icon("circle-notch", class = "fa-spin"), " ",
        if (ph$phase == "queued") "Waiting for a free worker…" else "Running…"))
    }
    if (!is.null(rv$run_error)) {
      div(class = "status-strip alert alert-danger",
          icon("circle-exclamation"), " ", rv$run_error)
    } else if (!is.null(rv$results)) {
      div(class = "status-strip alert alert-success",
          icon("circle-check"), " Analysis complete.")
    }
  })

  # ── Sidebar: currently loaded dataset (top of sidebar) ───────────────────
  output$current_dataset_ui <- renderUI({
    if (is.null(rv$imported_data)) {
      return(div(class = "current-dataset-box no-dataset",
        icon("circle-question"),
        span(class = "current-dataset-label", "No dataset loaded")
      ))
    }
    label <- if (is.null(rv$dataset_label)) "Custom dataset" else rv$dataset_label
    div(class = "current-dataset-box",
      icon("file-circle-check"),
      span(class = "current-dataset-label", label)
    )
  })

  # ── Sidebar: data load status (in Data section) ─────────────────────────
  output$data_load_status_ui <- renderUI({
    if (!is.null(rv$import_error)) {
      div(class = "status-strip alert alert-danger",
          icon("circle-exclamation"), " ", rv$import_error)
    } else if (!is.null(rv$imported_data)) {
      div(class = "status-strip alert alert-info",
          icon("circle-info"), " Data loaded.")
    }
  })

  # ── Sidebar: run status strip ────────────────────────────────────────────

  # ── Sidebar: download buttons ────────────────────────────────────────────
  output$download_ui <- renderUI({
    req(rv$results)
    tagList(
      hr(),
      p(class = "sidebar-section", icon("download"), " Download"),
      downloadButton("dl_posterior",  "Posterior probabilities", class = "btn-dl"),
      downloadButton("dl_comparison", "Bayesian vs match-counting", class = "btn-dl"),
      downloadButton("dl_convergence","Convergence summary",      class = "btn-dl")
    )
  })

  # Classification columns in both CSVs reflect the current threshold slider
  output$dl_posterior <- downloadHandler(
    filename = "posterior_probabilities.csv",
    content  = function(f) write.csv(posterior_results(), f, row.names = FALSE)
  )
  output$dl_comparison <- downloadHandler(
    filename = "bayesian_match_counting_comparison.csv",
    content  = function(f) write.csv(comparison_results(), f, row.names = FALSE)
  )
  output$dl_convergence <- downloadHandler(
    filename = "mcmc_convergence_summary.csv",
    content  = function(f) write.csv(rv$results$convergence, f, row.names = FALSE)
  )

  # ════════════════════════════════════════════════════════════════════════════
  #  Panel 1 — Explore Data
  # ════════════════════════════════════════════════════════════════════════════

  output$data_summary_ui <- renderUI({
    if (is.null(rv$imported_data)) {
      return(div(class = "welcome-box mb-3",
        p("Upload your genotyping data and marker metadata files using the sidebar to get started.")
      ))
    }

    d          <- rv$imported_data
    n_patients <- nrow(d$late_failures) / 2
    size_note  <- local({
      iter <- if (is.null(input$iter) || is.na(input$iter)) DEFAULT_MCMC$iter else input$iter
      problem <- run_size_problem(dataset_size()$patients, dataset_size()$additional, iter)
      if (!is.null(problem))
        div(class = "alert alert-warning mb-3", icon("triangle-exclamation"), " ", problem,
            " You can explore the data here; to classify it, split the dataset by site or use the MalReBay R package for larger studies (see About).")
    })
    n_markers  <- nrow(d$marker_info)
    n_sites    <- length(unique(d$late_failures$Site))

    tagList(div(class = "row mb-3",
      div(class = "col-sm-3",
        div(class = "stat-box",
          div(class = "stat-num teal", n_patients),
          div(class = "stat-label", "Patients")
        )
      ),
      div(class = "col-sm-3",
        div(class = "stat-box",
          div(class = "stat-num teal", n_sites),
          div(class = "stat-label", "Sites")
        )
      ),
      div(class = "col-sm-3",
        div(class = "stat-box",
          div(class = "stat-num teal", n_markers),
          div(class = "stat-label", "Markers")
        )
      ),
      div(class = "col-sm-3",
        div(class = "stat-box",
          div(class = "stat-num", style = "font-size:1.1rem; padding-top:0.35rem;",
              span(class = "badge badge-primary", d$data_type)),
          div(class = "stat-label", "Data type")
        )
      )
    ),
    size_note)
  })

  # Day 0 + recurrence rows, plus any extra Day 0 samples from the second sheet
  all_genotype_data <- reactive({
    req(rv$imported_data)
    dplyr::bind_rows(rv$imported_data$late_failures, rv$imported_data$additional)
  })

  # Site drop-downs of the Explore Data charts: "All sites" plus each site,
  # reset whenever another dataset is loaded
  observe({
    req(rv$imported_data)
    choices <- c(ALL_SITES, sort(unique(as.character(all_genotype_data()$Site))))
    for (id in c("loaded_data_site", "diversity_site", "allele_freq_site", "moi_site")) {
      updateSelectInput(session, id, choices = choices, selected = ALL_SITES)
    }
  })

  # Loaded Data: the recurrence (Day 0 + recurrence) sheet, or the
  # additional Day 0 sheet when the dataset has one
  has_additional_data <- reactive({
    additional <- rv$imported_data$additional
    !is.null(additional) && nrow(additional) > 0
  })

  output$loaded_data_set_ui <- renderUI({
    if (!has_additional_data()) return(NULL)
    radioButtons("loaded_data_set", "Data",
                 choices = c("Recurrence data" = "late_failures",
                             "Additional data" = "additional"))
  })

  output$loaded_data_table_ui <- renderUI({
    if (is.null(rv$imported_data)) {
      return(empty_panel_msg("table", "Upload or select a dataset to see the loaded data."))
    }
    DTOutput("loaded_data_table")
  })

  loaded_data_df <- reactive({
    req(rv$imported_data)
    # The switch's last value outlives the switch itself when the next
    # dataset has no additional data, so only honour it while it is shown
    use_additional <- has_additional_data() && identical(input$loaded_data_set, "additional")
    set <- if (use_additional) "additional" else "late_failures"
    df  <- rv$imported_data[[set]]
    req(!is.null(df))
    site <- input$loaded_data_site
    if (!is.null(site) && site %in% df$Site) df <- df[df$Site == site, , drop = FALSE]
    df
  })

  output$loaded_data_table <- renderDT({
    datatable(
      loaded_data_df(),
      options  = list(scrollX = TRUE, pageLength = 10),
      rownames = FALSE,
      class    = "table-sm table-striped table-bordered"
    )
  # Send the rows with the table rather than paging them via Ajax: datasets
  # are small, and a pending Ajax page request from the previous dataset's
  # table fails ("Invalid JSON response") once another dataset is loaded
  }, server = FALSE)

  # Pie charts are laid out <= 4 per row (plot_markers_diversity()'s max_cols
  # default), one per marker, like the allele distribution panels below
  DIVERSITY_MAX_COLS <- 4

  output$diversity_plot_ui <- renderUI({
    if (is.null(rv$imported_data)) {
      return(empty_panel_msg("chart-pie",
        "Upload or select a dataset to see allele / haplotype diversity."))
    }
    n_markers <- length(allele_dist_markers())
    req(n_markers > 0)
    n_row <- ceiling(n_markers / min(DIVERSITY_MAX_COLS, n_markers))
    plotOutput("diversity_plot", height = paste0(n_row * 280 + 40, "px"))
  })

  # list(all_sites = <figure>, by_site = list(<site> = <figure>)); built once
  # per dataset, so switching site only picks another figure
  diversity_plots <- reactive({
    # req() must come first: plot_markers_diversity() forces its data_type
    # argument before genotypedata, so a NULL imported_data reaches its
    # validation check (and throws a visible error) before genotypedata's
    # own req(rv$imported_data) ever gets a chance to run.
    req(rv$imported_data)
    plot_markers_diversity(
      genotypedata  = all_genotype_data(),
      data_type     = rv$imported_data$data_type,
      marker_info   = rv$imported_data$marker_info,
      output_folder = NULL,
      max_cols      = DIVERSITY_MAX_COLS
    )
  })

  output$diversity_plot <- renderPlot({
    plots <- diversity_plots()
    validate(need(!is.null(plots), "No diversity data available."))
    p <- select_site_plot(plots, input$diversity_site)
    validate(need(!is.null(p), "No diversity data available for this site."))
    print(p)
  })

  # Per-marker allele distribution from MalReBay::plot_allele_distribution():
  # histograms of fragment sizes for length-polymorphic markers, frequency
  # per haplotype for AmpSeq markers. One panel per marker, <= 3 per row.
  ALLELE_DIST_MAX_COLS <- 3

  # Markers with allele columns in the data that are also in the marker
  # metadata - the same set plot_allele_distribution() draws a panel for.
  # Cheap to compute, so the renderUI below can size the plot without
  # building it (the slow part stays in renderPlot, where the busy spinner
  # shows).
  allele_dist_markers <- reactive({
    req(rv$imported_data)
    marker_cols <- grep("(_allele_\\d+|_\\d+)$", colnames(all_genotype_data()), value = TRUE)
    intersect(unique(gsub("(_allele_\\d+|_\\d+)$", "", marker_cols)),
              rv$imported_data$marker_info$marker_id)
  })

  output$allele_freq_plot_ui <- renderUI({
    if (is.null(rv$imported_data)) {
      return(empty_panel_msg("chart-simple",
        "Upload or select a dataset to see allele frequency distribution."))
    }
    n_markers <- length(allele_dist_markers())
    req(n_markers > 0)
    n_row <- ceiling(n_markers / min(ALLELE_DIST_MAX_COLS, n_markers))
    plotOutput("allele_freq_plot", height = paste0(n_row * 300, "px"))
  })

  allele_freq_plots <- reactive({
    req(rv$imported_data)
    plot_allele_distribution(
      genotypedata  = all_genotype_data(),
      marker_info   = rv$imported_data$marker_info,
      output_folder = NULL,
      max_cols      = ALLELE_DIST_MAX_COLS
    )
  })

  output$allele_freq_plot <- renderPlot({
    plots <- allele_freq_plots()
    validate(need(!is.null(plots), "No allele frequency data available."))
    p <- select_site_plot(plots, input$allele_freq_site)
    validate(need(!is.null(p), "No allele frequency data available for this site."))
    print(p)
  })

  # plot_moi() returns a named list of one ggplot per site and has no pooled
  # version, so build that one by relabelling every sample's site as
  # "All sites" first.
  moi_plots <- reactive({
    data <- all_genotype_data()
    list(
      all_sites = plot_moi(dplyr::mutate(data, Site = ALL_SITES), output_folder = NULL)[[1]],
      by_site   = plot_moi(data, output_folder = NULL)
    )
  })

  output$moi_plot_ui <- renderUI({
    if (is.null(rv$imported_data)) {
      return(empty_panel_msg("layer-group",
        "Upload or select a dataset to see multiplicity of infection (MOI)."))
    }
    plotOutput("moi_plot", height = "400px")
  })

  output$moi_plot <- renderPlot({
    p <- select_site_plot(moi_plots(), input$moi_site)
    validate(need(!is.null(p), "No MOI data available."))
    p
  })

  # ════════════════════════════════════════════════════════════════════════════
  #  Panel 2 — Classification Results
  # ════════════════════════════════════════════════════════════════════════════

  # Elapsed time in words: "45 seconds", "2 min 14 s", "1 h 5 min"
  format_duration <- function(seconds) {
    seconds <- round(seconds)
    if (seconds < 60)   return(paste(seconds, if (seconds == 1) "second" else "seconds"))
    if (seconds < 3600) return(sprintf("%d min %d s", seconds %/% 60, seconds %% 60))
    sprintf("%d h %d min", seconds %/% 3600, (seconds %% 3600) %/% 60)
  }

  output$results_summary_ui <- renderUI({
    if (!is.null(rv$results) && !is.null(rv$run_info)) {
      info <- rv$run_info
      return(div(class = "run-info mb-3",
        icon("clock"), " ",
        sprintf("The MalReBay algorithm ran in %s (%d chains of %s iterations).",
                format_duration(info$seconds), info$n_chains,
                format(info$iter, big.mark = ","))))
    }
    if (is.null(rv$results)) {
      return(div(class = "welcome-box mb-3",
        p("Click", tags$strong("Run MalReBay"), "in the sidebar to run the Bayesian MCMC analysis.")
      ))
    }
    NULL
  })

  # Per-site breakdown: patients uploaded vs. recurrences actually classified
  # (these can differ if a patient's data couldn't be used, e.g. no
  # comparable loci), split into recrudescence / new infection counts.
  site_stats <- reactive({
    req(rv$results, rv$imported_data)

    patients_per_site <- rv$imported_data$late_failures %>%
      dplyr::mutate(Base.ID = trimws(gsub(" Day 0| recurrence", "", Sample.ID))) %>%
      dplyr::distinct(Site, Base.ID) %>%
      dplyr::count(Site, name = "Patients")

    recurrences_per_site <- posterior_results() %>%
      dplyr::group_by(Site) %>%
      dplyr::summarise(
        Recurrences    = dplyr::n(),
        Recrudescences = sum(Classification == "Recrudescence"),
        `New Infections` = sum(Classification == "New Infection"),
        .groups = "drop"
      )

    dplyr::full_join(patients_per_site, recurrences_per_site, by = "Site") %>%
      dplyr::mutate(Threshold = input$classification_threshold) %>%
      dplyr::arrange(Site)
  })

  output$site_stats_table_ui <- renderUI({
    if (is.null(rv$results)) {
      return(empty_panel_msg("table", "Run MalReBay to see per-site summary statistics."))
    }
    DTOutput("site_stats_table")
  })

  SITE_STATS_COLUMN_INFO <- c(
    Site             = "Study site.",
    Patients         = "Number of patients with genotyping data at this site.",
    Recurrences      = "Number of recurrences classified by MalReBay.",
    Recrudescences   = "Recurrences classified as recrudescence at the classification threshold.",
    `New Infections` = "Recurrences classified as new infection at the classification threshold.",
    Threshold        = "Classification threshold used: the value set in the sidebar."
  )

  output$site_stats_table <- renderDT({
    req(site_stats())
    datatable(
      site_stats(),
      container = header_with_info(names(site_stats()), SITE_STATS_COLUMN_INFO),
      options  = list(dom = "t", pageLength = 50),
      rownames = FALSE,
      class    = "table-sm table-striped table-bordered"
    )
  })

  # Site drop-downs of the results charts, reset after every run
  observe({
    req(rv$results)
    sites <- sort(unique(as.character(rv$results$posterior_probabilities$Site)))
    updateSelectInput(session, "prob_hist_site",
                      choices = c(ALL_SITES, EACH_SITE, sites), selected = ALL_SITES)
    updateSelectInput(session, "posterior_site",
                      choices = c(ALL_SITES, sites), selected = ALL_SITES)
    updateSelectInput(session, "comparison_site",
                      choices = c(ALL_SITES, sites), selected = ALL_SITES)
    conv_sites <- names(rv$results$mcmc_loglikelihoods)
    updateSelectInput(session, "convergence_site",
                      choices = conv_sites, selected = conv_sites[1])
  })

  output$prob_histogram_ui <- renderUI({
    if (is.null(rv$results)) {
      return(empty_panel_msg("chart-bar",
        "Run MalReBay to see the posterior probability distribution."))
    }
    # The faceted view gets 280px per row of site panels (the package lays
    # them out in a ceiling(sqrt(n_sites))-column grid)
    height <- 380
    if (identical(input$prob_hist_site, EACH_SITE)) {
      n_sites <- length(unique(rv$results$posterior_probabilities$Site))
      height  <- max(height, ceiling(n_sites / ceiling(sqrt(n_sites))) * 280 + 80)
    }
    div(style = "position:relative;",
      plotOutput("prob_histogram", height = paste0(height, "px"),
                 hover = hoverOpts(id = "prob_histogram_hover", delay = 50,
                                   delayType = "debounce", nullOutside = TRUE)),
      uiOutput("prob_histogram_tooltip")
    )
  })

  # MalReBay::plot_probability_histogram() with the slider's threshold: the
  # pooled histogram, the faceted per-site one, or (for a single site) the
  # pooled histogram of that site's patients only
  output$prob_histogram <- renderPlot({
    req(rv$results, input$classification_threshold)
    pp <- rv$results$posterior_probabilities
    req(nrow(pp) > 0)
    site <- input$prob_hist_site

    single_site <- !is.null(site) && site %in% pp$Site
    if (single_site) pp <- pp[pp$Site == site, , drop = FALSE]

    plots <- without_printing(plot_probability_histogram(
      list(posterior_probabilities = pp),
      output_folder  = NULL,
      verbose        = FALSE,
      prob_threshold = input$classification_threshold
    ))
    req(plots)

    if (single_site) {
      plots$all_sites + labs(title = paste("Posterior Probability Distribution -", site))
    } else if (identical(site, EACH_SITE)) {
      plots$by_site
    } else {
      plots$all_sites
    }
  })

  # Hover tooltip: probabilities near 0.5 (a fixed 0.4-0.6 band, independent
  # of the adjustable classification threshold) carry the most uncertainty,
  # regardless of which side of the threshold they land on.
  output$prob_histogram_tooltip <- renderUI({
    req(input$prob_histogram_hover)
    x <- input$prob_histogram_hover$x
    req(!is.na(x), x >= 0, x <= 1)

    text <- if (x >= 0.4 && x <= 0.6) "Less confidence in classification" else "Strong confidence in classification"
    css  <- input$prob_histogram_hover$coords_css

    div(style = paste0(
          "position:absolute; left:", css$x + 12, "px; top:", css$y - 10, "px; ",
          "background:#2C3E50; color:#fff; padding:4px 8px; border-radius:4px; ",
          "font-size:0.78rem; font-weight:600; pointer-events:none; z-index:20; white-space:nowrap;"
        ),
        text)
  })

  output$posterior_table_ui <- renderUI({
    if (is.null(rv$results)) {
      return(empty_panel_msg("table", "Run MalReBay to see patient-level results."))
    }
    DTOutput("posterior_table")
  })

  # White (new infection) -> blue -> red (recrudescence) diverging
  # "Outcome" scale.
  OUTCOME_WHITE <- "#F7F7F7"
  OUTCOME_BLUE  <- "#67A9CF"
  OUTCOME_RED   <- "#D6604D"
  prob_gradient_colors <- grDevices::colorRampPalette(
    c(OUTCOME_WHITE, "#ade8f4", OUTCOME_BLUE, "#F4A582", OUTCOME_RED)
  )(100)
  prob_gradient_breaks <- seq(0, 1, length.out = 101)[2:100]

  classify_prob <- function(p) {
    ifelse(p >= input$classification_threshold, "Recrudescence", "New Infection")
  }

  # posterior_probabilities with Classification re-derived from the slider
  # (MalReBay() only applied the threshold that was set at run time) and
  # moved next to Probability; the N_Markers_* columns follow.
  posterior_results <- reactive({
    req(rv$results)
    rv$results$posterior_probabilities %>%
      dplyr::mutate(Classification = classify_prob(Probability),
                    Threshold      = input$classification_threshold) %>%
      dplyr::relocate(Threshold, Classification, .after = Probability)
  })

  # Hover text for the Patient Results column headers
  POSTERIOR_COLUMN_INFO <- c(
    Site                 = "Study site of the patient.",
    Sample.ID            = "Patient identifier (shared by the Day 0 and recurrence samples).",
    Probability          = "Posterior probability that the recurrence is a recrudescence (0 to 1).",
    Classification       = "Recrudescence if the probability is at or above the classification threshold set in the sidebar, otherwise new infection.",
    Threshold            = "Classification threshold used: the value set in the sidebar.",
    N_Markers_Day0       = "Number of markers with genotyping data in the Day 0 sample.",
    N_Markers_Recurrence = "Number of markers with genotyping data in the recurrence sample.",
    N_Markers_Compared   = "Number of markers typed in both samples, which are the ones compared between Day 0 and recurrence."
  )

  # A DT header row with a hover tooltip (title attribute) on each column
  header_with_info <- function(columns, info) {
    tags$table(tags$thead(tags$tr(lapply(columns, function(col) {
      tags$th(title = if (col %in% names(info)) info[[col]], col)
    }))))
  }

  output$posterior_table <- renderDT({
    pp <- posterior_results()
    site <- input$posterior_site
    if (!is.null(site) && site != ALL_SITES) pp <- pp[pp$Site == site, , drop = FALSE]
    pp$Probability <- round(pp$Probability, 4)

    datatable(
      pp,
      container = header_with_info(names(pp), POSTERIOR_COLUMN_INFO),
      options  = list(scrollX = TRUE, pageLength = 15, order = list(list(2, "desc"))),
      rownames = FALSE,
      class    = "table-sm table-striped table-bordered"
    ) |>
      formatStyle(
        "Classification",
        color      = styleEqual(c("Recrudescence", "New Infection"), c(OUTCOME_RED, OUTCOME_BLUE)),
        fontWeight = "bold"
      ) |>
      formatStyle(
        "Probability",
        backgroundColor = styleInterval(prob_gradient_breaks, prob_gradient_colors),
        fontWeight       = "600"
      )
  })

  # MalReBay's comparison table: one row per patient with the per-marker
  # match-counting calls (R/NI/IND/ERR), combined msp1/msp2 family calls for
  # MSP panels, the overall loose/strict rules (e.g. Match_counting_5of7 /
  # Match_counting_7of7, or 2of3 / 3of3 for MSP panels) as
  # "Recrudescence"/"New infection"/NA, then MalReBay_probability and
  # MalReBay_classification. Classification is re-derived from the slider.
  comparison_results <- reactive({
    req(rv$results, rv$results$comparison)
    rv$results$comparison %>%
      dplyr::mutate(MalReBay_classification = classify_prob(MalReBay_probability),
                    MalReBay_threshold      = input$classification_threshold) %>%
      dplyr::relocate(MalReBay_threshold, .after = MalReBay_probability)
  })

  # Display version: overall match-counting rules recoded to the same
  # R/NI/IND shorthand as the per-marker calls (NA = too few markers scored
  # to apply the rule) and readable headers for those rules. MalReBay's own
  # columns keep their MalReBay_ prefix, to tell them apart from the
  # match-counting calls.
  match_counting_summary <- reactive({
    df <- comparison_results()

    rule_cols  <- grep("^Match_counting_", colnames(df), value = TRUE)
    recode_who <- function(x) dplyr::case_when(x == "Recrudescence" ~ "R",
                                               x == "New infection" ~ "NI",
                                               TRUE ~ "IND")
    df[rule_cols] <- lapply(df[rule_cols], recode_who)

    df$MalReBay_probability <- round(df$MalReBay_probability, 4)
    # "Match_counting_5of7" -> "Match counting 5/7"
    colnames(df) <- sub("^Match_counting_(\\d+)of(\\d+)$", "Match counting \\1/\\2", colnames(df))
    df
  })

  output$comparison_table_ui <- renderUI({
    if (is.null(rv$results)) {
      return(empty_panel_msg("table",
        "Run MalReBay to compare the MalReBay and match-counting classifications."))
    }
    DTOutput("comparison_table")
  })

  # Hover text for the comparison table's column headers. Marker and rule
  # columns depend on the panel, so their text is built from the names
  # ("Match counting 5/7": at least 5 of the 7 markers must match).
  comparison_column_info <- function(call_cols) {
    rule_cols   <- grep("^Match counting ", call_cols, value = TRUE)
    marker_cols <- setdiff(call_cols, rule_cols)

    marker_info <- setNames(paste0(
      "Match-counting result at marker ", marker_cols, ": R if the Day 0 and ",
      "recurrence samples share an allele, NI if they do not, IND if the marker ",
      "was not typed in both samples."), marker_cols)

    rule_info <- vapply(rule_cols, function(col) {
      n <- as.integer(strsplit(sub("^Match counting ", "", col), "/")[[1]])
      needed <- if (n[1] == n[2]) paste0("all ", n[2], " markers match")
                else paste0("at least ", n[1], " of the ", n[2], " markers match")
      paste0("Overall match-counting classification: R (recrudescence) if ", needed,
             ", otherwise NI (new infection). IND if fewer than ", n[1],
             " markers could be compared.")
    }, character(1))

    c(Site      = "Study site of the patient.",
      Sample.ID = "Patient identifier (shared by the Day 0 and recurrence samples).",
      marker_info, rule_info,
      MalReBay_probability    = "MalReBay posterior probability that the recurrence is a recrudescence (0 to 1).",
      MalReBay_classification = "MalReBay classification: recrudescence if the probability is at or above the classification threshold set in the sidebar, otherwise new infection.",
      MalReBay_threshold      = "Classification threshold used for the MalReBay classification: the value set in the sidebar.")
  }

  output$comparison_table <- renderDT({
    df        <- match_counting_summary()
    site      <- input$comparison_site
    if (!is.null(site) && site != ALL_SITES) df <- df[df$Site == site, , drop = FALSE]
    call_cols <- setdiff(colnames(df), c("Site", "Sample.ID", "MalReBay_probability",
                                         "MalReBay_classification", "MalReBay_threshold"))
    prob_idx  <- which(colnames(df) == "MalReBay_probability") - 1  # DT columns are 0-based

    dt <- datatable(
      df,
      container = header_with_info(names(df), comparison_column_info(call_cols)),
      options  = list(scrollX = TRUE, pageLength = 15, order = list(list(prob_idx, "desc"))),
      rownames = FALSE,
      class    = "table-sm table-striped table-bordered"
    ) |>
      formatStyle(
        "MalReBay_classification",
        color      = styleEqual(c("Recrudescence", "New Infection"), c(OUTCOME_RED, OUTCOME_BLUE)),
        fontWeight = "bold"
      ) |>
      formatStyle(
        "MalReBay_probability",
        backgroundColor = styleInterval(prob_gradient_breaks, prob_gradient_colors),
        fontWeight       = "600"
      )

    # Per-marker calls as plain centred text; the overall match-counting
    # classifications in the same red / blue as recrudescence / new infection
    # elsewhere, with indeterminate and error calls in grey
    rule_cols   <- grep("^Match counting ", call_cols, value = TRUE)
    marker_cols <- setdiff(call_cols, rule_cols)
    if (length(marker_cols) > 0) {
      dt <- dt |> formatStyle(marker_cols, textAlign = "center")
    }
    for (col in rule_cols) {
      dt <- dt |> formatStyle(
        col,
        backgroundColor = styleEqual(c("R", "NI", "IND", "ERR"),
                                     c(OUTCOME_RED, OUTCOME_BLUE, "#F0F0F0", "#F0F0F0")),
        color            = styleEqual(c("R", "NI", "IND", "ERR"),
                                     c("#FFFFFF", "#FFFFFF", "#95A5A6", "#95A5A6")),
        fontWeight       = "700",
        textAlign        = "center"
      )
    }

    dt
  })

  # MalReBay::plot_comparison_heatmap(): match-counting rules vs MalReBay's
  # posterior probability, one heatmap per site plus one with all sites
  # (patients grouped by site). Uses the probability itself, not the
  # threshold, so it doesn't change with the slider.
  heatmap_plots <- reactive({
    req(rv$results, rv$imported_data)
    without_printing(plot_comparison_heatmap(
      rv$results,
      marker_info   = rv$imported_data$marker_info,
      output_folder = NULL,
      verbose       = FALSE
    ))
  })

  output$heatmap_plot_ui <- renderUI({
    if (is.null(rv$results)) {
      return(empty_panel_msg("table-cells",
        "Run MalReBay to see the match-counting vs MalReBay heatmap."))
    }
    # Tiles shrink with patient count; grow the plot so rows stay readable
    comp    <- rv$results$comparison
    site    <- input$comparison_site
    n_rows  <- if (!is.null(site) && site %in% comp$Site) sum(comp$Site == site) else nrow(comp)
    plotOutput("heatmap_plot", height = paste0(max(450, n_rows * 10 + 150), "px"))
  })

  output$heatmap_plot <- renderPlot({
    plots <- heatmap_plots()
    p <- select_site_plot(plots, input$comparison_site)
    validate(need(!is.null(p), "No comparison data available."))
    p
  })

  output$convergence_table_ui <- renderUI({
    if (is.null(rv$results)) {
      return(empty_panel_msg("gauge-high", "Run MalReBay to see MCMC convergence diagnostics."))
    }
    if (is.null(rv$results$convergence)) {
      return(empty_panel_msg("triangle-exclamation",
        "No convergence diagnostics are available for this run."))
    }
    tagList(
      p(class = "table-legend",
        "The table shows several criteria for checking that the MCMC chains have
         converged, i.e. that the results are reliable. Each criterion has its own
         threshold, and it can be acceptable for one or two of them to narrowly
         fail, for example a slightly low effective sample size. In that case,
         check the diagnostic plots below: well-mixed, overlapping chains indicate
         convergence. If several criteria fail clearly, consider re-running with
         more iterations. Hover over a column name for its meaning and threshold."),
      DTOutput("convergence_table"),
      p(class = "table-legend table-legend-below",
        span(class = "legend-swatch", style = paste0("background:", PASS_COLOR)), " criterion met ",
        span(class = "legend-swatch", style = paste0("background:", FAIL_COLOR)), " criterion not met")
    )
  })

  # Convergence criteria, as used by MalReBay's own diagnostics: R-hat close
  # to 1 (Gelman-Rubin < 1.1, rank-normalised < 1.01) and an effective sample
  # size above 400. All are computed on the log-likelihood of the chains.
  CONVERGENCE_CRITERIA <- list(
    Gelman_Rubin_Rhat     = function(x) x < 1.1,
    Rank_Rhat             = function(x) x < 1.01,
    ESS_Bulk              = function(x) x > 400,
    ESS_Tail              = function(x) x > 400,
    Effective_Sample_Size = function(x) x > 400
  )
  CONVERGENCE_COLUMN_INFO <- c(
    Site                  = "Study site. The MCMC is run separately for each site.",
    Gelman_Rubin_Rhat     = "Classical Gelman-Rubin R-hat: compares the variation between and within the chains. Passes below 1.1.",
    Rank_Rhat             = "Rank-normalised R-hat (Vehtari et al. 2021), a stricter and more robust version of R-hat. Passes below 1.01.",
    ESS_Bulk              = "Bulk effective sample size: the number of independent draws the chains are worth for estimating the centre of the distribution. Passes above 400.",
    ESS_Tail              = "Tail effective sample size: as bulk ESS, but for the tails of the distribution. Passes above 400.",
    Effective_Sample_Size = "Classical effective sample size (coda), summed over the chains. Passes above 400."
  )
  PASS_COLOR <- "#D4EDDA"
  FAIL_COLOR <- "#F8D7DA"

  output$convergence_table <- renderDT({
    req(rv$results, rv$results$convergence)
    conv <- rv$results$convergence
    metrics <- intersect(names(CONVERGENCE_CRITERIA), names(conv))

    # One hidden pass/fail column per metric drives that metric's cell colour
    # (missing values stay uncoloured)
    status_cols <- paste0(metrics, "_status")
    for (i in seq_along(metrics)) {
      ok <- CONVERGENCE_CRITERIA[[metrics[i]]](conv[[metrics[i]]])
      conv[[status_cols[i]]] <- ifelse(is.na(ok), "", ifelse(ok, "pass", "fail"))
    }
    shown <- setdiff(names(conv), status_cols)

    dt <- datatable(
      conv,
      container = header_with_info(names(conv), CONVERGENCE_COLUMN_INFO),
      options   = list(dom = "t", pageLength = 50,
                       columnDefs = list(list(visible = FALSE,
                                              targets = match(status_cols, names(conv)) - 1))),
      rownames  = FALSE,
      class     = "table-sm table-striped table-bordered"
    ) |>
      formatRound(columns = intersect(c("Gelman_Rubin_Rhat", "Rank_Rhat"), shown), digits = 4) |>
      formatRound(columns = intersect(c("ESS_Bulk", "ESS_Tail", "Effective_Sample_Size"), shown), digits = 0)

    for (i in seq_along(metrics)) {
      dt <- formatStyle(dt, metrics[i], valueColumns = status_cols[i],
                        backgroundColor = styleEqual(c("pass", "fail"), c(PASS_COLOR, FAIL_COLOR)))
    }
    dt
  })

  output$convergence_plots_ui <- renderUI({
    if (is.null(rv$results)) {
      return(empty_panel_msg("images", "Run MalReBay to see convergence diagnostic plots."))
    }
    if (length(rv$results$mcmc_loglikelihoods) == 0) {
      return(empty_panel_msg("triangle-exclamation",
        "No MCMC chains are available for convergence plots."))
    }
    plotOutput("convergence_plot", height = "950px")
  })

  # The four diagnostic panels of MalReBay::plot_likelihood_diagnostics() for
  # the selected site, rebuilt from the run's per-chain log-likelihoods
  # (combine_plots = TRUE arranges them 2x2; no Gelman-Rubin with one chain)
  output$convergence_plot <- renderPlot({
    req(rv$results, input$convergence_site)
    chains <- rv$results$mcmc_loglikelihoods[[input$convergence_site]]
    req(chains)
    diag <- plot_likelihood_diagnostics(
      all_chains_loglikelihood = chains,
      site_name                = input$convergence_site,
      save_plot                = FALSE,
      verbose                  = FALSE,
      combine_plots            = TRUE
    )
    validate(need(!is.null(diag), "No valid chains for this site."))
    diag$combined
  })

  # ════════════════════════════════════════════════════════════════════════════
  #  Panel 3 — Drug Efficacy
  # ════════════════════════════════════════════════════════════════════════════
  # PCR-corrected efficacy needs the total per-protocol enrolled population per
  # site (including patients who never recurred), which classification results
  # alone don't contain - only recurrences are genotyped/imported. The
  # enrollment inputs below default to the observed recurrence count (the
  # minimum possible enrollment) and must be raised to the true site total for
  # an accurate estimate.

  output$efficacy_summary_ui <- renderUI({
    if (is.null(rv$results)) {
      return(div(class = "welcome-box mb-3",
        p("Click", tags$strong("Run MalReBay"), "in the sidebar to compute
          PCR-corrected drug efficacy from the classification results.")
      ))
    }
    NULL
  })

  output$efficacy_enrollment_inputs_ui <- renderUI({
    if (is.null(rv$results)) {
      return(empty_panel_msg("user-group", "Run MalReBay to enter per-site enrollment."))
    }
    stats <- site_stats()
    req(stats)

    # Default enrollment = total patients recorded anywhere in the dataset
    # for that site: recurrence-patients (late_failures) plus any Day-0-only
    # background samples (the "additional" sheet). Used only until the user
    # enters the real enrolled number - Enrolled == Recurrences alone would
    # force PCR-corrected efficacy to exactly 0% for every site (every
    # recurrence is by definition either a recrudescence or reinfection, so
    # the numerator is always zero), a guaranteed-misleading first look.
    additional <- rv$imported_data$additional
    additional_per_site <- if (!is.null(additional) && nrow(additional) > 0) {
      table(additional$Site)
    } else {
      integer(0)
    }

    rows <- lapply(seq_len(nrow(stats)), function(i) {
      site     <- stats$Site[i]
      recur    <- stats$Recurrences[i]
      if (is.na(recur)) recur <- 0
      patients <- stats$Patients[i]
      if (is.na(patients)) patients <- 0

      extra <- if (as.character(site) %in% names(additional_per_site)) {
        additional_per_site[[as.character(site)]]
      } else {
        0
      }
      default_enrolled <- max(patients + extra, recur)

      div(class = "row mb-2 align-items-center",
        div(class = "col-sm-6", tags$strong(site),
          span(class = "text-muted", style = "font-size:0.8rem;",
               paste0(" (", recur, " recurrence(s); ", default_enrolled,
                      " patients in dataset)"))),
        div(class = "col-sm-6",
          numericInput(paste0("enrolled_", make.names(site)),
                       label = NULL, value = default_enrolled, min = recur, step = 1)
        )
      )
    })
    do.call(tagList, rows)
  })

  # Site-level table combining classification counts with the user-entered
  # enrollment numbers and the resulting PCR-corrected efficacy:
  #   efficacy % = (enrolled - recrudescences) / (enrolled - reinfections) * 100
  # NA when enrollment hasn't been entered yet or is below the observed
  # recurrence count (an impossible enrollment given the data).
  efficacy_data <- reactive({
    req(rv$results)
    stats <- site_stats()
    req(stats)

    stats$Recurrences     <- ifelse(is.na(stats$Recurrences),     0, stats$Recurrences)
    stats$Recrudescences  <- ifelse(is.na(stats$Recrudescences),  0, stats$Recrudescences)
    stats$`New Infections` <- ifelse(is.na(stats$`New Infections`), 0, stats$`New Infections`)

    stats$Enrolled <- vapply(stats$Site, function(site) {
      val <- input[[paste0("enrolled_", make.names(site))]]
      if (is.null(val) || is.na(val)) NA_real_ else as.numeric(val)
    }, numeric(1))

    # Enrolled == Recurrences is the bare minimum the input allows (nobody
    # enrolled beyond those who recurred), which forces Efficacy to exactly
    # 0% - every recurrence is by definition either a recrudescence or
    # reinfection, so the numerator collapses to zero. That's a real,
    # if degenerate, answer for that scenario; the input's default value is
    # deliberately set above this minimum so it isn't the first thing shown.
    valid <- !is.na(stats$Enrolled) & stats$Enrolled >= stats$Recurrences
    stats$Efficacy <- NA_real_
    stats$Efficacy[valid] <- 100 *
      (stats$Enrolled[valid] - stats$`New Infections`[valid] - stats$Recrudescences[valid]) /
      pmax(stats$Enrolled[valid] - stats$`New Infections`[valid], 1)

    stats
  })

  output$efficacy_table_ui <- renderUI({
    if (is.null(rv$results)) {
      return(empty_panel_msg("pills", "Run MalReBay to see drug efficacy estimates."))
    }
    DTOutput("efficacy_table")
  })

  output$efficacy_table <- renderDT({
    df <- efficacy_data()
    req(df)

    display <- df %>%
      dplyr::transmute(
        Site,
        Recurrences,
        `Recrudescence (true failure)` = Recrudescences,
        `Reinfection`                  = `New Infections`,
        `Enrolled`                     = Enrolled,
        `PCR-Corrected Efficacy (%)`   = round(Efficacy, 1)
      )

    datatable(
      display,
      options  = list(dom = "t", pageLength = 50),
      rownames = FALSE,
      class    = "table-sm table-striped table-bordered"
    ) |>
      formatStyle(
        "PCR-Corrected Efficacy (%)",
        color      = styleInterval(90, c("#D6604D", "#1E8449")),
        fontWeight = "700"
      )
  })

  output$efficacy_plot_ui <- renderUI({
    if (is.null(rv$results)) {
      return(empty_panel_msg("chart-column", "Run MalReBay to see efficacy by site."))
    }
    plotOutput("efficacy_plot", height = "360px")
  })

  output$efficacy_plot <- renderPlot({
    df <- efficacy_data()
    req(df)
    df <- df[!is.na(df$Efficacy), , drop = FALSE]
    validate(need(nrow(df) > 0,
      "Enter a valid enrolled population (>= observed recurrences) for at least one site."))

    ggplot(df, aes(x = Site, y = Efficacy, fill = Efficacy >= 90)) +
      geom_col(width = 0.55) +
      geom_hline(yintercept = 90, linetype = "dashed", color = "#7F8C8D") +
      geom_text(aes(label = paste0(round(Efficacy, 1), "%")), vjust = -0.6, size = 4.2) +
      scale_fill_manual(values = c(`TRUE` = "#2C3E50", `FALSE` = "#D6604D"), guide = "none") +
      labs(x = NULL, y = "PCR-corrected efficacy (%)") +
      ylim(0, 105) +
      theme_minimal(base_size = 13) +
      theme(panel.grid.major.x = element_blank())
  })
}
