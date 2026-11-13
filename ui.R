# ── "How It Works" page building blocks ──────────────────────────────────────

# Click-to-open information box: a small (i) button showing `...` in a
# popover, so the page carries one-line summaries and the detail on demand
hiw_info <- function(title, ...) {
  popover(
    tags$button(type = "button", class = "hiw-info-btn",
                `aria-label` = paste("More about:", title), icon("circle-info")),
    ..., title = title
  )
}

# A control with its label and info box side by side. The info button sits
# outside the input's own <label> (inputs get label = NULL), so clicking it
# doesn't also focus or open the input.
hiw_field <- function(label, info, input) {
  div(class = "hiw-field",
    div(class = "hiw-field-label", span(label), if (!is.null(info)) hiw_info(label, info)),
    input
  )
}

# One step of the "reasoning behind the result" flow: a one-line summary on
# the card; the whole card opens the full explanation when clicked
hiw_flow_step <- function(number, title, summary, details, final = FALSE) {
  popover(
    div(class = paste("hiw-flow-step", if (final) "hiw-flow-step-final"),
        role = "button", tabindex = "0",
      div(class = "hiw-flow-head", span(class = "hiw-step-badge", number), title),
      p(summary),
      span(class = "hiw-more", icon("circle-info"), " Learn more")
    ),
    details, title = title
  )
}

# Slider whose raw numbers mean nothing to a non-technical reader: the
# numeric labels are hidden in CSS and replaced by "Low" / "High" ends.
hiw_qualitative_slider <- function(id, label, value, min = 0.05, max = 0.95, step = 0.05,
                                   info = NULL) {
  hiw_field(label, info,
    div(class = "hiw-qual-slider",
      sliderInput(id, NULL, min = min, max = max, value = value, step = step, ticks = FALSE),
      div(class = "hiw-qual-ends", span("Low"), span("High"))
    )
  )
}

# Most alleles/haplotypes detectable per sample in the toy example, i.e. a
# multiplicity of infection (MOI) of up to 3. Kept small so the pairs table
# and chart labels stay readable (3 x 3 = 9 pairs); the real model has no
# such cap.
HIW_MAX_MOI <- 3

# Up to HIW_MAX_MOI alleles/haplotypes picked at one timepoint
hiw_allele_picker <- function(id, label, choices, selected) {
  hiw_field(label, NULL,
    selectizeInput(id, NULL, choices = choices, selected = selected, multiple = TRUE,
                   options = list(maxItems = HIW_MAX_MOI,
                                  placeholder = paste("Choose 1 to", HIW_MAX_MOI))))
}

# A numbered panel with its controls placed right next to its chart
hiw_panel <- function(number, title, controls, plot, caption, info = NULL) {
  div(class = "content-card",
    div(class = "content-card-header", span(class = "hiw-step-badge", number), " ", title,
        if (!is.null(info)) hiw_info(title, info)),
    div(class = "hiw-card-body",
      div(class = "row",
        div(class = "col-md-4 hiw-controls", controls),
        div(class = "col-md-8", plot, p(class = "hiw-caption", caption))
      )
    )
  )
}

# Right-hand result panel: strength of evidence -> starting assumption ->
# probability of recrudescence. Sticky, so it stays in view next to the
# charts that drive it.
hiw_result_panel <- function(prefix) {
  div(class = "content-card hiw-result-card",
    div(class = "content-card-header", "Evidence and result"),
    div(class = "hiw-card-body",
      uiOutput(paste0(prefix, "_evidence_ui")),
      hr(),
      hiw_field("Starting assumption: chance of recrudescence",
        "The chance of recrudescence before looking at the genetic data (the
         prior). MalReBay uses a neutral 0.5; move it to see how much the
         starting point matters compared with the evidence.",
        sliderInput(paste0(prefix, "_prior"), NULL,
                    min = 0.05, max = 0.95, value = 0.5, step = 0.05, ticks = FALSE)),
      hr(),
      # Gauge and verdict are separate outputs around a fixed slider, so the
      # slider isn't re-created (and reset) every time the probability changes
      uiOutput(paste0(prefix, "_prob_ui")),
      hiw_field("Classification threshold",
        "The probability at or above which the case is classified as a
         recrudescence (dashed line on the bar above). MalReBay uses 0.5 by
         default; a higher threshold asks for stronger evidence before calling
         a treatment failure.",
        sliderInput(paste0(prefix, "_threshold"), NULL,
                    min = 0.5, max = 0.95, value = 0.5, step = 0.05, ticks = FALSE)),
      uiOutput(paste0(prefix, "_verdict_ui"))
    )
  )
}

# Comparison with the match-counting method, toggles next to the table
hiw_compare_panel <- function(prefix) {
  match_toggle <- function(n, selected) {
    hiw_field(paste("Marker", n), NULL,
      radioButtons(paste0(prefix, "_match_", n), NULL, inline = TRUE,
                   choices = c("Match" = "match", "Mismatch" = "mismatch"),
                   selected = selected))
  }

  div(class = "content-card",
    div(class = "content-card-header", "Compare: match-counting ",
      hiw_info("Match-counting",
        "The standard method counts how many markers match exactly between
         Day 0 and recurrence, and calls a recrudescence if enough of them
         match (all of them for the strict rule, most of them for the loose
         rule). Unlike MalReBay, it ignores how common the alleles are and how
         close two non-identical alleles are.")),
    div(class = "hiw-card-body",
      div(class = "row",
        div(class = "col-md-5 hiw-controls",
          p(class = "hiw-help", "Marker 1 is the example above. Set two more
             markers to complete the panel."),
          match_toggle(2, "match"),
          match_toggle(3, "mismatch")
        ),
        div(class = "col-md-7", uiOutput(paste0(prefix, "_match_counting_ui")))
      )
    )
  )
}

# ── About page building blocks ──────────────────────────────────────────────

# Partner logos for the About page, from www/logos/. A missing file is
# skipped rather than shown as a broken image.
PARTNER_LOGOS <- list(
  list(file = "SwissTPH_Logo_Screen_Indigo_rgb.png", alt = "Swiss Tropical and Public Health Institute",
       href = "https://www.swisstph.ch"),
  list(file = "Gates_Foundation_Logo.svg",           alt = "Gates Foundation",
       href = "https://www.gatesfoundation.org")
)

partner_logos <- function(which = seq_along(PARTNER_LOGOS)) {
  logos <- lapply(PARTNER_LOGOS[which], function(logo) {
    if (!file.exists(file.path("www", "logos", logo$file))) return(NULL)
    tags$a(href = logo$href, target = "_blank", class = "about-logo",
           tags$img(src = paste0("logos/", logo$file), alt = logo$alt))
  })
  logos <- Filter(Negate(is.null), logos)
  if (length(logos) > 0) div(class = "about-logos", logos)
}

# WHO overview page on monitoring antimalarial drug efficacy (Welcome page)
WHO_DRUG_EFFICACY_URL <- paste0("https://www.who.int/teams/global-malaria-programme/",
  "case-management/drug-efficacy-and-resistance/monitoring-antimalarial-drug-efficacy")

# MalReBay package resources linked from the About page
MALREBAY_GITHUB <- "https://github.com/SwissTPH/MalReBay"
MALREBAY_SITE   <- "https://swisstph.github.io/MalReBay"

about_link <- function(icon_name, label, href, description) {
  tags$li(
    icon(icon_name, class = "about-link-icon"),
    div(tags$a(href = href, target = "_blank", label), span(description))
  )
}

# ── Analysis page building blocks ────────────────────────────────────────────

# A chart card: short description, an optional site drop-down (choices are
# filled in by the server once data / results are available), then the chart
plot_card <- function(icon_name, title, description, site_input_id, output_id,
                      collapsible = FALSE, open = TRUE, extra_controls = NULL) {
  body <- div(style = "padding:1rem;",
    div(class = "card-intro-row",
      div(class = "card-intro", description),
      div(class = "card-controls",
        if (!is.null(site_input_id))
          div(class = "card-site-select",
            selectInput(site_input_id, "Choose site:", choices = ALL_SITES, width = "220px")),
        extra_controls)
    ),
    uiOutput(output_id)
  )
  header_content <- tagList(icon(icon_name), " ", title)

  if (!collapsible) {
    return(div(class = "content-card",
      div(class = "content-card-header", header_content), body))
  }
  collapsible_card(header_content, body, open = open)
}

# A content card whose header opens and closes its body: a native <details>
# element. www/collapse.js tells Shiny when a card is opened so its outputs
# render or resume.
collapsible_card <- function(header_content, body, open = TRUE) {
  tags$details(class = "content-card collapsible-card", open = if (open) NA,
    tags$summary(class = "content-card-header",
      header_content, icon("chevron-down", class = "collapse-chevron")),
    body)
}

ui <- navbarPage(
  title  = tagList(
    tags$img(src = "MalReBay_logo_1.svg", alt = "MalReBay logo", class = "navbar-logo"),
    div(class = "navbar-title-block",
      span(class = "navbar-title", "MalReBay"),
      span(class = "navbar-subtitle", "Bayesian Classification of Malaria Recurrence")
    )
  ),
  # Browser tab title (otherwise derived from the title tags above, subtitle included)
  windowTitle = "MalReBay",
  id     = "main_nav",
  theme  = bs_theme(bootswatch = "flatly"),
  header = tagList(
    useShinyjs(),
    # Spinning circle centred over any plot/table while it's being computed
    useBusyIndicators(spinners = TRUE, pulse = FALSE),
    busyIndicatorOptions(spinner_type  = "ring",
                         spinner_color = "#2C3E50",
                         spinner_size  = "40px",
                         spinner_delay = "0.2s"),
    tags$head(
      # ?v=<modification time>: browsers reload the stylesheet after it changes
      # instead of serving a cached copy
      tags$link(rel = "stylesheet",
                href = paste0("custom.css?v=", as.integer(file.mtime(file.path("www", "custom.css"))))),
      tags$script(src = "collapse.js"),
      # Reports the content height to the parent page when embedded in an iframe
      tags$script(src = "iframe-height.js"),
      tags$script(HTML("
        $(document).ready(function() {
          $(document).on('click', '#load_example_lp_short, #load_example_lp, #load_example_ampseq', function() {
            $('#load_example_lp_short, #load_example_lp, #load_example_ampseq').removeClass('btn-example-selected');
            $(this).addClass('btn-example-selected');
          });
        });
      "))
    )
  ),

  # ── Welcome / Landing page ────────────────────────────────────────────────
  tabPanel(
    title = "Welcome",

    div(class = "lp-body",

      # ── Row 1: MalReBay method ────────────────────────────────────────────
      div(class = "lp-card lp-card-wide",
        h4(class = "lp-card-title", "About MalReBay"),
        p("MalReBay estimates, for each patient, the ",
          tags$strong("posterior probability that a recurrence is a recrudescence"),
          " using a fully Bayesian statistical model."),
        div(class = "lp-feature-grid",
          div(class = "lp-feature",
            div(
              tags$strong("Handles polyclonal infections"),
              p("Models patients carrying multiple parasite clones at Day 0
                or at recurrence.")
            )
          ),
          div(class = "lp-feature",
            div(
              tags$strong("Imputes missing genotypes"),
              p("Missing allele calls are treated as latent variables and inferred
                jointly with the classification probability.")
            )
          ),
          div(class = "lp-feature",
            div(
              tags$strong("Supports multiple marker types"),
              p("Length-polymorphic microsatellites, MSP/GLURP size variants, and
                amplicon sequencing (AmpSeq) haplotypes.")
            )
          ),
          div(class = "lp-feature",
            div(
              tags$strong("Bayesian Sampling and Reliability Checks"),
              p("Runs a Bayesian analysis and checks whether the results are stable
                and reliable.")
            )
          )
        ),
        div(class = "lp-cta",
          actionButton("go_to_home", "Start Analysis",
                       class = "btn btn-primary lp-cta-btn")
        )
      ),

      # ── Row 2: TES + WHO ──────────────────────────────────────────────────
      div(class = "lp-row",

        div(class = "lp-card",
          h4(class = "lp-card-title", "Therapeutic Efficacy Studies"),
          p("Therapeutic Efficacy Studies (TES) are the gold-standard field trials used to
            monitor the clinical and parasitological response to antimalarial drugs.
            Patients with confirmed ", tags$em("P. falciparum"), " malaria are treated and
            then followed for 28–42 days. Any parasitaemia detected during follow-up is
            classified as either:"),
          tags$ul(
            tags$li(tags$strong("Recrudescence"), " — the original parasite clone persists,
              indicating treatment failure;"),
            tags$li(tags$strong("Reinfection"), " — a new parasite clone acquired after
              treatment, not a sign of drug failure.")
          ),
          p("Distinguishing recrudescence from reinfection requires molecular
            genotyping, because both events are clinically indistinguishable.")
        ),

        div(class = "lp-card",
          h4(class = "lp-card-title", "WHO Guidelines"),
          p("The World Health Organization recommends conducting TES at sentinel sites
            every two years to detect emerging antimalarial resistance. The ",
            tags$a(href = WHO_DRUG_EFFICACY_URL, target = "_blank",
                   "WHO guidance on monitoring antimalarial drug efficacy"),
            " includes:"),
          tags$ul(
            tags$li("Enrolment of patients with uncomplicated ",
              tags$em("P. falciparum"), " monoinfection;"),
            tags$li("28-day or 42-day follow-up with standardised outcome classification;"),
            # One HTML string: separate tags$em() children would be rendered
            # on their own lines, leaving a visible space before each comma
            tags$li(HTML("Molecular genotyping of <em>pfmsp1</em>, <em>pfmsp2</em> and
              a panel of microsatellite markers (replacing <em>pfglurp</em>) to
              distinguish recrudescence from new infection, with amplicon
              sequencing (AmpSeq) encouraged where available;")),
            tags$li("A site treatment failure rate ≥ 10 % triggers a policy review.")
          ),
          p("MalReBay offers a Bayesian approach to the analysis of the resulting
            genotyping data: it compares the Day 0 and day of recurrence samples
            and classifies each recurrence as a recrudescence or a new infection,
            as an alternative to the standard allele match-counting method.")
        )
      ),

      # ── Row 3: Authors ────────────────────────────────────────────────────
      div(class = "lp-card lp-card-wide lp-authors-card",
        h4(class = "lp-card-title", icon("users", class = "lp-icon-dark"), " Authors"),
        div(class = "lp-authors",
          div(class = "lp-author",
            icon("user-circle", class = "lp-author-icon"),
            div(
              tags$strong("Monica Golumbeanu"),
              p("Swiss Tropical and Public Health Institute (Swiss TPH), Basel, Switzerland")
            )
          ),
          div(class = "lp-author",
            icon("user-circle", class = "lp-author-icon"),
            div(
              tags$strong("Veronica Adhiambo Ochieng"),
              p("Swiss Tropical and Public Health Institute (Swiss TPH), Basel, Switzerland"),
              p("African Institute for Mathematical Sciences (AIMS), Kigali, Rwanda")
            )
          )
        ),
        partner_logos(which = 1),  # Swiss TPH
        p(class = "lp-github",
          icon("github"), " Source code: ",
          tags$a(href = "https://github.com/SwissTPH/MalReBay", target = "_blank",
                 "github.com/SwissTPH/MalReBay")),
        p(class = "lp-funding",
          "The development of MalReBay and the associated dashboard is funded
           by the Gates Foundation.")
      )
    )
  ),

  # ── How It Works (self-contained toy example, no Stan/MalReBay calls) ─────
  tabPanel(
    title = "How It Works",

    div(class = "hiw-page",

      h1(class = "lp-title", "How MalReBay Classifies a Recurrence"),
      p(class = "hiw-intro",
        "When malaria occurs again after treatment, is it because the original
         infection was not fully cleared (recrudescence), or because the person
         has acquired a new infection? MalReBay identifies this by comparing the
         alleles detected in the parasite population of the Day 0 sample with
         those detected in the recurrence sample. It then calculates a
         likelihood score to assess whether the data are more consistent with a
         recrudescence or a new infection."),

      # ── The reasoning, in plain terms ──────────────────────────────────────
      div(class = "content-card",
        div(class = "content-card-header", "The reasoning in four steps ",
          span(class = "hiw-header-hint", "(click a step for details)")),
        div(class = "hiw-card-body",
          div(class = "hiw-flow",
            hiw_flow_step("1", "How common?",
              "A shared rare allele is strong evidence.",
              "At each genetic marker a parasite carries one variant, called an
               allele. If Day 0 and recurrence share an allele that is rare
               locally, a new infection is unlikely to carry it by chance: strong
               evidence of recrudescence. A shared common allele is weak
               evidence."),
            div(class = "hiw-flow-join", "+"),
            hiw_flow_step("2", "How similar?",
              "Closer alleles point to the same parasite.",
              "Identical alleles point to the same parasite. Lab measurements are
               not perfectly precise, so slightly different alleles may still
               come from the same parasite, but the larger the difference, the
               weaker the link."),
            div(class = "hiw-flow-join", icon("arrow-right")),
            hiw_flow_step("3", "Strength of evidence",
              "Steps 1 + 2 combined into one number.",
              "The likelihood ratio: how many times more likely the observed
               alleles are under recrudescence than under a new infection. Above
               1 supports recrudescence; below 1 supports a new infection."),
            div(class = "hiw-flow-join", icon("arrow-right")),
            hiw_flow_step("4", "Probability",
              "From a neutral 0.5 to a final verdict.",
              tagList(
                p("The evidence updates a neutral starting probability (0.5) into a
                   probability of recrudescence between 0 and 1. At or above
                   the threshold (0.5 by default) the case is called a
                   recrudescence."),
                p(class = "mb-0", "In a real analysis this is repeated across all
                   markers, and patients carrying several parasite strains at
                   once are accounted for.")
              ),
              final = TRUE)
          )
        )
      ),

      # ── Interactive example ────────────────────────────────────────────────
      h3(class = "hiw-section-title", "Try it yourself"),
      p(class = "hiw-intro",
        "Let's consider a single genetic marker with ten possible alleles, and
         one patient whose Day 0 and recurrence samples were both genotyped at
         this marker. Choose which alleles were detected in each sample, how
         diverse the local parasite population is and how precise the
         measurement is, and see how the strength of the evidence and the
         probability of recrudescence on the right change."),
      p(class = "hiw-intro",
        paste0("For simplicity, each sample can carry at most ", HIW_MAX_MOI,
               " alleles here, i.e. a multiplicity of infection (MOI) of up to ",
               HIW_MAX_MOI, ". MalReBay itself handles infections with more clones.")),
      radioButtons("hiw_data_type", "Choose the type of genetic data",
        choices = c("Length-polymorphic markers (microsatellites, MSP1/MSP2/GLURP)" = "lp",
                    "Amplicon sequencing (microhaplotypes)"                        = "ampseq"),
        selected = "lp", width = "100%"),

      # ── Length-polymorphic: alleles are sizes; similarity = size difference
      conditionalPanel("input.hiw_data_type == 'lp'",
        div(class = "row",
          div(class = "col-lg-8",
            hiw_panel("1", "How common are the alleles?",
              controls = tagList(
                hiw_allele_picker("hiw_lp_d0", "Day 0 alleles",
                                  as.character(1:10), "4"),
                hiw_allele_picker("hiw_lp_dx", "Recurrence alleles",
                                  as.character(1:10), c("4", "7")),
                hiw_qualitative_slider("hiw_lp_evenness", "Local genetic diversity",
                  value = 0.55,
                  info = "Move the slider to the right to make the alleles more
                          equally common (high diversity): a shared allele is then
                          rarer and more telling. Move it to the left to let a few
                          alleles dominate (low diversity): sharing one of them
                          says little.")
              ),
              plot    = plotOutput("hiw_lp_freq_plot", height = "300px"),
              caption = "Bars: how common each allele size is locally. Triangles:
                         this patient's alleles."
            ),
            hiw_panel("2", "How similar are the alleles?",
              controls = tagList(
                hiw_qualitative_slider("hiw_lp_dparam", "Measurement precision",
                  value = 0.6,
                  info = "Move the slider to the right for more precise sizing:
                          the bars concentrate on a difference of 0, so only
                          identical or nearly identical sizes support
                          recrudescence. Move it to the left to accept larger size
                          differences as measurement error.")
              ),
              plot    = plotOutput("hiw_lp_kernel_plot", height = "240px"),
              caption = "Each bar is the likelihood of that size difference if it
                         is a recrudescence (same parasite). Identical sizes
                         (difference 0) are the most likely. Each further step is
                         less likely than the one before by the same factor, which
                         is set by the measurement precision. The bars add up to 1.
                         Red bars and markers: the differences between this
                         patient's Day 0 and recurrence alleles.",
              info = tagList(
                p("If the recurrence is a recrudescence, the same parasite is
                   measured twice, so its Day 0 and recurrence allele sizes should
                   be identical. Sizing is not perfect, though, so small
                   differences can occur."),
                p("The chart shows how likely each size difference is under
                   recrudescence. A difference of 0 is the most likely, and each
                   extra step is less likely than the one before by a fixed
                   factor, so the bars decrease steadily. They add up to 1."),
                p("Measurement precision sets how quickly they decrease. For each
                   pair of a Day 0 and a recurrence allele, the bar at their size
                   difference is that pair's L(R), the likelihood if
                   recrudescence used in step 3.")
              )
            )
          ),
          div(class = "col-lg-4", hiw_result_panel("hiw_lp"))
        ),
        hiw_compare_panel("hiw_lp")
      ),

      # ── Amplicon sequencing: alleles are sequences; similarity = exact match
      conditionalPanel("input.hiw_data_type == 'ampseq'",
        div(class = "row",
          div(class = "col-lg-8",
            hiw_panel("1", "How common are the haplotypes?",
              controls = tagList(
                hiw_allele_picker("hiw_hap_d0", "Day 0 haplotypes",
                                  paste0("H", 1:10), "H4"),
                hiw_allele_picker("hiw_hap_dx", "Recurrence haplotypes",
                                  paste0("H", 1:10), c("H4", "H7")),
                hiw_qualitative_slider("hiw_hap_evenness", "Local genetic diversity",
                  value = 0.55,
                  info = "Move the slider to the right to make the haplotypes
                          more equally common (high diversity): a shared haplotype
                          is then rarer and more telling. Move it to the left to
                          let a few haplotypes dominate (low diversity): sharing
                          one of them says little.")
              ),
              plot    = plotOutput("hiw_hap_freq_plot", height = "300px"),
              caption = "Bars: how common each haplotype is locally. Triangles:
                         this patient's haplotypes."
            ),
            hiw_panel("2", "Do the haplotypes match?",
              controls = tagList(
                hiw_qualitative_slider("hiw_hap_qq", "Sequencing error allowance",
                  value = 0.08, min = 0.01, max = 0.4, step = 0.01,
                  info = "Sequenced haplotypes either match exactly or not. This
                          is the small chance that the same parasite is read as a
                          different haplotype. Move the slider to the right to
                          allow more error: a mismatch then counts less against
                          recrudescence. Move it to the left to make a mismatch
                          more decisive.")
              ),
              plot    = plotOutput("hiw_hap_error_plot", height = "240px"),
              caption = "Likelihood of each outcome if it is a recrudescence (same
                         parasite). The same haplotype is read again with
                         likelihood 1 − error allowance. The remaining chance is
                         shared equally among the other 9 haplotypes, so a specific
                         different haplotype gets error allowance ÷ 9.
                         Red bars and markers: the comparisons between this
                         patient's Day 0 and recurrence haplotypes.",
              info = tagList(
                p("If the recurrence is a recrudescence, the same parasite is
                   sequenced twice, so the haplotype should be read identically:
                   likelihood 1 − error allowance."),
                p("With the remaining chance (the error allowance), it is misread
                   as one of the other 9 haplotypes, each equally likely: error
                   allowance ÷ 9 for a specific different haplotype."),
                p("For each pair of a Day 0 and a recurrence haplotype, the
                   matching bar is that pair's L(R), the likelihood if
                   recrudescence used in step 3.")
              )
            )
          ),
          div(class = "col-lg-4", hiw_result_panel("hiw_hap"))
        ),
        hiw_compare_panel("hiw_hap")
      )
    )
  ),

  # ── Analysis page ────────────────────────────────────────────────────────
  tabPanel(
    title = "Analysis",

    # Custom flex layout replaces sidebarLayout so we can animate sidebar width
    div(class = "dashboard-layout",

      # ── Sidebar ──────────────────────────────────────────────────────────
      div(id = "sidebar-col",

        # Collapse / expand toggle button
        div(class = "sidebar-toggle-bar",
          actionButton("sidebar_toggle", NULL,
                       icon  = icon("chevron-left"),
                       class = "btn-sidebar-toggle",
                       title = "Collapse sidebar")
        ),

        # All sidebar content (hidden when collapsed)
        div(id = "sidebar-content",

          # ── Currently loaded dataset ─────────────────────────────────────
          uiOutput("current_dataset_ui"),

          # ── LOAD DATA ──────────────────────────────────────────────────
          div(class = "sidebar-main-section", "Load data"),

          # ── Example datasets (collapsible, closed by default) ───────────
          tags$a(
            class            = "sidebar-subsection mcmc-toggle d-flex justify-content-between align-items-center",
            href             = "#examples-collapse",
            `data-bs-toggle` = "collapse",
            `aria-expanded`  = "false",
            `aria-controls`  = "examples-collapse",
            style            = "text-decoration:none; cursor:pointer;",
            span(icon("database"), " Example datasets"),
            icon("chevron-down", class = "mcmc-chevron")
          ),

          tags$div(
            id    = "examples-collapse",
            class = "collapse",

            actionButton("load_example_lp_short",
                         span(icon("vial"), " Length Polymorphic - 1 site"),
                         class = "btn-example"),
            actionButton("load_example_lp",
                         span(icon("vial"), " Length Polymorphic - 2 sites"),
                         class = "btn-example"),
            actionButton("load_example_ampseq",
                         span(icon("dna"), " Amplicon Sequencing - 1 site"),
                         class = "btn-example")
          ),

          uiOutput("data_load_status_ui"),

          fileInput("genotype_file",
            label = span("Genotyping data ",
              tooltip(icon("circle-info", class = "info-icon"),
                "Excel (.xlsx) file with one Day 0 row and one recurrence row per patient.
                 Allele columns are named <marker_id>_1, <marker_id>_2, etc.
                 Numeric values for length-polymorphic markers; haplotype strings for AmpSeq.",
                placement = "right")),
            placeholder = ".xlsx",
            accept      = c(".xlsx", ".xls"),
            buttonLabel = "Browse"),

          # Blank input templates bundled with MalReBay, one per data type.
          # data-bs-display="static" keeps the menu inside the sidebar, whose
          # overflow:hidden would otherwise clip a Popper-positioned menu.
          div(class = "file-helper dropdown",
            tags$a(class            = "file-helper-link dropdown-toggle",
                   href             = "#",
                   role             = "button",
                   `data-bs-toggle` = "dropdown",
                   `data-bs-display` = "static",
                   `aria-expanded`  = "false",
                   icon("file-arrow-down"), " Input data templates"),
            tags$ul(class = "dropdown-menu file-helper-menu",
              tags$li(static_download_link("MalReBay_template_microsatellites.xlsx",
                                           "Microsatellite panel", class = "dropdown-item")),
              tags$li(static_download_link("MalReBay_template_msp_glurp.xlsx",
                                           "MSP1 / MSP2 / GLURP", class = "dropdown-item")),
              tags$li(static_download_link("MalReBay_template_ampseq.xlsx",
                                           "Amplicon sequencing", class = "dropdown-item"))
            )
          ),

          fileInput("marker_file",
            label = span("Marker metadata ",
              tooltip(icon("circle-info", class = "info-icon"),
                "Excel (.xlsx) file with one row per marker. Required columns:
                 marker_id, marker_type, repeatlength, binning_method.",
                placement = "right")),
            placeholder = ".xlsx",
            accept      = c(".xlsx", ".xls"),
            buttonLabel = "Browse"),

          div(class = "file-helper",
            static_download_link("makers_details.xlsx",
                                 span(icon("file-arrow-down"), " Marker metadata example"),
                                 class = "file-helper-link")
          ),

          # ── CLASSIFICATION ─────────────────────────────────────────────
          div(class = "sidebar-main-section", "Classification"),

          sliderInput("classification_threshold",
            label = span("Recrudescence threshold ",
              tooltip(icon("circle-info", class = "info-icon"),
                "Posterior probability of recrudescence at or above this
                 value is classified as Recrudescence, below it as New
                 Infection. Default: 0.5.",
                placement = "right")),
            min = 0.5, max = 1, value = 0.5, step = 0.01),

          # ── MCMC settings (collapsible, closed by default) ──────────────
          tags$a(
            class            = "sidebar-subsection mcmc-toggle d-flex justify-content-between align-items-center",
            href             = "#mcmc-collapse",
            `data-bs-toggle` = "collapse",
            `aria-expanded`  = "false",
            `aria-controls`  = "mcmc-collapse",
            style            = "text-decoration:none; cursor:pointer;",
            span(icon("sliders"), " MCMC Settings"),
            icon("chevron-down", class = "mcmc-chevron")
          ),

          tags$div(
            id    = "mcmc-collapse",
            class = "collapse",

            numericInput("n_chains",
              label = span("Chains ",
                tooltip(icon("circle-info", class = "info-icon"),
                  "Number of independent Stan MCMC chains.
                   More chains improve convergence diagnostics but increase run time.
                   Default and maximum: 3.",
                  placement = "right")),
              value = DEFAULT_MCMC$n_chains, min = 1, max = MAX_MCMC_CHAINS, step = 1),

            numericInput("iter",
              label = span("Total iterations ",
                tooltip(icon("circle-info", class = "info-icon"),
                  paste0("Total number of iterations per chain (warm-up + sampling combined).
                   Split into warm-up and sampling according to the burn-in fraction below.
                   Default: 2000. Maximum: ", format(MAX_ITERATIONS, big.mark = ","), "."),
                  placement = "right")),
              value = DEFAULT_MCMC$iter, min = 500, max = MAX_ITERATIONS, step = 500),

            numericInput("burn_in_frac",
              label = span("Burn-in fraction ",
                tooltip(icon("circle-info", class = "info-icon"),
                  "Proportion of total iterations used as Stan warm-up (adaptation).
                   E.g. 0.5 uses the first 50 % of iterations as warm-up.
                   Default: 0.5.",
                  placement = "right")),
              value = DEFAULT_MCMC$burn_in_frac, min = 0.05, max = 0.90, step = 0.05),

            numericInput("adapt_delta",
              label = span("Adapt delta ",
                tooltip(icon("circle-info", class = "info-icon"),
                  "Stan's target average proposal acceptance probability.
                   Higher values (e.g. 0.95-0.99) reduce divergent transitions
                   at the cost of slower sampling. Default: 0.85.",
                  placement = "right")),
              value = DEFAULT_MCMC$adapt_delta, min = 0.5, max = 0.999, step = 0.01),

            numericInput("random_seed",
              label = span("Random seed ",
                tooltip(icon("circle-info", class = "info-icon"),
                  "Seed for the MCMC sampler, for reproducible results.
                   Default: 42.",
                  placement = "right")),
              value = DEFAULT_MCMC$random_seed, min = 1, step = 1)
          ),

          actionButton("run_btn",
                       label = span(icon("play"), " Run MalReBay"),
                       class = "btn-primary"),
          p(class = "run-limit-note",
            sprintf("Dashboard limit: %d patients, %s iterations. Split by site or use the MalReBay R package for larger studies.",
                    MAX_RECURRENT_PATIENTS, format(MAX_ITERATIONS, big.mark = ","))),

          uiOutput("download_ui"),
          uiOutput("run_status_ui")
        )
      ),

      # ── Main panel ────────────────────────────────────────────────────────
      div(id = "main-col",

        # Panel toggle bar
        div(class = "panel-toggle-bar",
          actionButton("btn_data_panel",
                       span(icon("chart-pie"), " Explore Data"),
                       class = "btn-panel-toggle active-panel"),
          actionButton("btn_results_panel",
                       span(icon("circle-check"), " Classification Results"),
                       class = "btn-panel-toggle"),
          actionButton("btn_efficacy_panel",
                       span(icon("pills"), " Drug Efficacy"),
                       class = "btn-panel-toggle")
        ),

        # ── Panel 1: Explore Data ─────────────────────────────────────────
        div(id = "data_exploration_panel",

          uiOutput("data_summary_ui"),

          p(class = "panel-intro",
            "A descriptive analysis of the genotyping data:
             how diverse each marker is, how its alleles are distributed, and
             how many parasite clones the samples carry. Each chart pools Day 0
             and recurrence samples (plus any additional Day 0 samples) and can
             be shown for all sites together or for a single site."),

          plot_card("table", "Loaded Data",
            p("Genotyping data are imported with one row per sample and one
               column for each detected allele at each marker (numbered _1, _2,
               ... until the maximum observed MOI for the respective marker).
               Cells are left empty when fewer alleles are detected in a sample
               for a given marker. The recurrence data contain the paired Day 0
               and recurrence samples. The additional data, if provided, contain
               extra Day 0 samples."),
            "loaded_data_site", "loaded_data_table_ui", collapsible = TRUE, open = FALSE,
            extra_controls = uiOutput("loaded_data_set_ui")),

          plot_card("chart-pie", "Allele / Haplotype Diversity",
            tagList(
              p("One pie chart per marker showing how often each allele (binned
                 fragment size for length-polymorphic markers) or haplotype
                 (amplicon sequencing) is observed. ", tags$em("n"), " is the
                 number of infections typed at that marker."),
              p("Many small slices mean high diversity: a match between Day 0
                 and recurrence is then unlikely to happen by chance, so the
                 marker is informative. One dominant slice means low diversity
                 and weaker evidence from that marker.")
            ),
            "diversity_site", "diversity_plot_ui", collapsible = TRUE),

          plot_card("chart-simple", "Allele Frequency Distribution",
            tagList(
              p("The distribution of raw allele values for each marker, as a
                 percentage of that marker's calls. For length-polymorphic
                 markers this is a histogram of fragment sizes (bp) and for
                 amplicon sequencing it is the frequency of each haplotype.")
            ),
            "allele_freq_site", "allele_freq_plot_ui", collapsible = TRUE),

          plot_card("layer-group", "Multiplicity of Infection (MOI)",
            tagList(
              p("The number of distinct alleles detected per marker in each
                 sample, an indicator of how many parasite clones an infection
                 carries. Violins show the distribution, the inner box the
                 median and interquartile range, points the individual samples,
                 and the number above each marker the mean MOI.")
            ),
            "moi_site", "moi_plot_ui", collapsible = TRUE)
        ),

        # ── Panel 2: Classification Results ──────────────────────────────
        div(id = "classification_panel", style = "display:none;",

          uiOutput("results_summary_ui"),

          collapsible_card(
            tagList(icon("table"), " Per-Site Summary",
              tags$span(id = "spinner_site_stats", class = "tab-spinner",
                icon("circle-notch", class = "fa-spin"))),
            div(style = "padding:1rem;",
              div(class = "card-intro-row",
                div(class = "card-intro",
                  p("An overview of the classification results for each site."))),
              uiOutput("site_stats_table_ui"))
          ),

          collapsible_card(
            tagList(icon("chart-bar"), " Probability of Recrudescence Distribution",
              tags$span(id = "spinner_prob_histogram", class = "tab-spinner",
                icon("circle-notch", class = "fa-spin"))),
            div(style = "padding:1rem;",
              div(class = "card-intro-row",
                div(class = "card-intro",
                  p("Number of patients by posterior probability of recrudescence,
                     coloured by classification at the current threshold (dotted
                     line, set in the sidebar). Most patients near 0 or 1 means
                     the markers were informative; values near the threshold are
                     the least certain. The plot shows at a glance how clearly
                     the patients separate into recrudescences and new
                     infections, and how many classifications depend on the
                     choice of threshold.")),
                div(class = "card-site-select",
                  selectInput("prob_hist_site", "Choose site:", choices = ALL_SITES, width = "220px"))
              ),
              uiOutput("prob_histogram_ui"))
          ),

          collapsible_card(
            tagList(icon("list-check"), " Classification Results"),
            div(class = "inner-tabs", style = "padding:1rem;",
            tabsetPanel(
              type = "tabs",

              tabPanel(
                span("Patient Results",
                  tags$span(id = "spinner_patient_results", class = "tab-spinner",
                    icon("circle-notch", class = "fa-spin"))),
                br(),
                div(class = "card-site-select mb-3",
                  selectInput("posterior_site", "Choose site:", choices = ALL_SITES, width = "220px")),
                div(class = "card-intro mb-2",
                  p("Each row represents one patient and shows the MalReBay
                     posterior probability of recrudescence, the resulting
                     classification at the current threshold, and the number of
                     markers available at Day 0, at recurrence and in both samples.
                     Hover over a column name for its meaning.")),
                uiOutput("posterior_table_ui")
              ),

              tabPanel(
                span("MalReBay vs Match-Counting",
                  tags$span(id = "spinner_comparison", class = "tab-spinner",
                    icon("circle-notch", class = "fa-spin"))),
                br(),
                # One site choice for both the table and the heatmap below, so
                # they always show the same patients
                div(class = "card-site-select mb-3",
                  selectInput("comparison_site", "Choose site:", choices = ALL_SITES, width = "220px")),
                plot_card("table-list", "Match-Counting vs MalReBay Results Table",
                  tagList(
                    p("Each row represents one patient and compares the MalReBay
                       result with the match-counting approach. It shows the match
                       result for each marker, the overall match-counting result
                       (e.g. whether 2 out of 3 or all 3 markers match), and the
                       corresponding MalReBay result."),
                    # One HTML string: separate tags$strong() children would
                    # leave a visible space before each colon
                    p(HTML("<strong>R</strong>: recrudescence (alleles match) &middot;
                       <strong>NI</strong>: new infection (no match) &middot;
                       <strong>IND</strong>: indeterminate (marker not typed in both
                       samples; for a rule, too few markers typed) &middot;
                       <strong>ERR</strong>: no unique Day 0 / recurrence pair."))),
                  NULL, "comparison_table_ui", collapsible = TRUE),
                br(),
                plot_card("table-cells", "Match-Counting vs MalReBay Heatmap",
                  p("One row per patient. The two match-counting rules are shown
                     as recrudescence (red) or new infection (pale blue); MalReBay
                     shows the posterior probability on the same colour scale.
                     Grey: too few markers scored to apply the rule."),
                  NULL, "heatmap_plot_ui", collapsible = TRUE)
              ),

              tabPanel(
                span("Convergence Diagnostics",
                  tags$span(id = "spinner_convergence", class = "tab-spinner",
                    icon("circle-notch", class = "fa-spin"))),
                br(),
                uiOutput("convergence_table_ui"),
                br(),
                div(class = "card-intro-row",
                  div(class = "card-intro",
                    p("Trace plot, Gelman-Rubin shrink factor, log-posterior
                       distribution and autocorrelation of the MCMC chains for the
                       selected site. Well-mixed, overlapping chains and a shrink
                       factor close to 1 indicate convergence.")),
                  div(class = "card-site-select",
                    selectInput("convergence_site", "Choose site:", choices = NULL, width = "220px"))
                ),
                uiOutput("convergence_plots_ui")
              )
            )
            )
          )
        ),

        # ── Panel 3: Drug Efficacy ────────────────────────────────────────
        div(id = "efficacy_panel", style = "display:none;",

          uiOutput("efficacy_summary_ui"),

          div(class = "content-card",
            div(class = "content-card-header", icon("user-group"), " Per-Protocol Enrollment",
              tooltip(icon("circle-info", class = "info-icon"),
                "Total number of patients enrolled and followed up at each site
                 (including those who never recurred). Needed to convert
                 recrudescence/reinfection counts into a PCR-corrected efficacy
                 rate. Defaults to the total number of patients recorded in
                 the dataset for that site - correct it to the true enrolled
                 number for an accurate estimate.",
                placement = "right")
            ),
            div(style = "padding:1rem;",
              p(class = "text-muted", style = "font-size:0.85rem;",
                "Default values are taken from the uploaded data. Adapt as
                 needed with the actual enrolled numbers for each site."),
              uiOutput("efficacy_enrollment_inputs_ui")
            )
          ),

          div(class = "content-card",
            div(class = "content-card-header", icon("table"), " PCR-Corrected Drug Efficacy"),
            div(style = "padding:1rem;", uiOutput("efficacy_table_ui"))
          ),

          div(class = "content-card",
            div(class = "content-card-header", icon("chart-column"), " Efficacy by Site"),
            div(style = "padding:1rem;", uiOutput("efficacy_plot_ui"))
          )
        )
      )
    )
  ),

  # ── About ─────────────────────────────────────────────────────────────────
  tabPanel(
    title = "About",

    div(class = "about-page",

      h1(class = "lp-title", "About MalReBay"),
      p(class = "about-lead",
        "MalReBay supports malaria researchers and national malaria programmes
         in distinguishing recrudescence of the original infection from new
         infection in patients with recurrent parasitaemia after antimalarial
         treatment. This classification underpins the calculation of
         PCR-corrected therapeutic efficacy, and MalReBay aims to make these
         estimates more robust and accurate."),

      div(class = "about-section",
        h4("What the dashboard can do"),
        p("The dashboard runs a complete MalReBay analysis in the browser, without
           installing R or writing code:"),
        tags$ul(class = "about-capabilities",
          tags$li(HTML("<strong>Load data</strong> from the Excel templates for
            microsatellites, MSP1/MSP2/GLURP or amplicon sequencing, or try one of the
            example datasets.")),
          tags$li(HTML("<strong>Explore the data</strong> before classification: the
            imported table, marker diversity, allele size distributions and the
            multiplicity of infection, for all sites together or per site.")),
          tags$li(HTML("<strong>Classify recurrences</strong> with the Bayesian model,
            site by site, with an adjustable classification threshold, a comparison with
            the match-counting method and MCMC convergence diagnostics.")),
          tags$li(HTML("<strong>Estimate PCR-corrected drug efficacy</strong> per site
            from the classification results and the number of patients enrolled.")),
          tags$li(HTML("<strong>Download</strong> the posterior probabilities, the
            match-counting comparison and the convergence summary as CSV files."))
        ),
        # One HTML string: separate tags$strong() children would leave a
        # visible space before each comma
        p(HTML(sprintf(
          "To keep the shared server responsive, an analysis in the dashboard is limited to
           <strong>%d patients with a recurrence</strong>, <strong>%s additional Day 0
           samples</strong>, <strong>%d MCMC chains</strong> and <strong>%s iterations per
           chain</strong>. Larger studies can be split by site and analysed separately, or
           analysed with the MalReBay R package below, which has no such limits.",
          MAX_RECURRENT_PATIENTS, format(MAX_ADDITIONAL_SAMPLES, big.mark = ","),
          MAX_MCMC_CHAINS, format(MAX_ITERATIONS, big.mark = ","))))
      ),

      div(class = "about-section",
        h4("The MalReBay R package"),
        p("This dashboard is built on MalReBay, an open-source R package. Using
           the package directly provides full control over the analysis and
           supports reproducible, documented workflows."),
        tags$ul(class = "about-links",
          about_link("house", "Package website", MALREBAY_SITE,
                     "Overview, installation instructions and the function reference."),
          about_link("book-open", "Vignette", paste0(MALREBAY_SITE, "/articles/MalReBay.html"),
                     "A complete walkthrough of an analysis, from data import to the
                      interpretation of results."),
          about_link("file-code", "Analysis notebooks",
                     paste0(MALREBAY_GITHUB, "/tree/main/inst/notebooks"),
                     "Focused notebooks for the analysis of microsatellite,
                      MSP1/MSP2/GLURP and amplicon sequencing data."),
          about_link("file-arrow-down", "Input data templates",
                     paste0(MALREBAY_GITHUB, "/tree/main/inst/templates"),
                     "Excel and CSV templates for each data type."),
          about_link("github", "Source code", MALREBAY_GITHUB,
                     "The package repository on GitHub, for issue reports and
                      contributions.")
        )
      ),

      div(class = "about-section",
        h4("Background"),
        p("Therapeutic efficacy studies follow patients for several weeks after
           treatment to assess the efficacy of antimalarial medicines. When
           parasites reappear during follow-up, only a recrudescence is counted
           as a treatment failure. The distinction is made by comparing the
           alleles detected in the parasite population before and after
           treatment, and the resulting efficacy estimates inform national
           treatment policies.")
      ),

      div(class = "about-section",
        h4("Approach"),
        p("MalReBay implements a Bayesian approach to this classification,
           which complements the established allele match-counting method. For
           each patient, it estimates the posterior probability that the
           recurrence is a recrudescence. The model accounts for the local
           frequency of each allele, the similarity between alleles, infections
           with multiple parasite clones and missing genotyping data. It
           supports microsatellite, MSP1/MSP2/GLURP and amplicon sequencing
           data."),
        p("This dashboard provides access to the method without running the
           package in R or requiring technical expertise.
           The ",
          actionLink("about_to_hiw", "How It Works"), " page illustrates the
           underlying reasoning with an interactive example.")
      ),

      div(class = "about-section",
        h4("Who we are"),
        p("MalReBay was developed at Swiss TPH in Basel, Switzerland."),
        div(class = "about-people",
          div(class = "about-person",
            tags$strong("Monica Golumbeanu"),
            span("Swiss TPH, Basel, Switzerland")),
          div(class = "about-person",
            tags$strong("Veronica Adhiambo Ochieng"),
            span("Swiss TPH, Basel, Switzerland; African Institute for
                  Mathematical Sciences (AIMS), Kigali, Rwanda"))
        ),
        partner_logos(which = 1)  # Swiss TPH
      ),

      div(class = "about-section",
        h4("Funding"),
        p("The development of MalReBay is funded by the Gates Foundation."),
        partner_logos(which = 2)  # Gates Foundation
      )
    )
  )
)

# navbarPage() emits a bare .navbar (no .navbar-expand-*), which bslib's
# Bootstrap 3 compat layer only lays out horizontally at >= 992px; below that
# the tabs stack into a tall vertical list (no hamburger, since the navbar
# isn't collapsible). That's what happens when the app is embedded in a
# narrow iframe on the mirror site, so force the horizontal layout at every
# width - custom.css shrinks it to fit narrow frames instead.
ui <- htmltools::tagQuery(ui)$find("nav.navbar")$addClass("navbar-expand")$allTags()
