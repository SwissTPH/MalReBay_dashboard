<p align="center">
  <img src="www/MalReBay_logo_1.svg" alt="MalReBay logo" width="140">
</p>

<h1 align="center">MalReBay Dashboard</h1>

<p align="center"><em>Bayesian Classification of Malaria Recurrence</em></p>

<p align="center">
  <a href="https://swisstph.shinyapps.io/malrebay-dashboard/"><strong>Open the dashboard</strong></a>
</p>

MalReBay supports malaria researchers and national malaria programmes in
distinguishing recrudescence of the original infection from new infection in
patients with recurrent parasitaemia after antimalarial treatment. This
classification underpins the calculation of PCR-corrected therapeutic
efficacy, and MalReBay aims to make these estimates more robust and accurate.

## What the dashboard can do

The dashboard runs a complete MalReBay analysis in the browser, without
installing R or writing code:

- **Load data** from the Excel templates for microsatellites, MSP1/MSP2/GLURP
  or amplicon sequencing, or try one of the example datasets.
- **Explore the data** before classification: the imported table, marker
  diversity, allele size distributions and the multiplicity of infection, for
  all sites together or per site.
- **Classify recurrences** with the Bayesian model, site by site, with an
  adjustable classification threshold, a comparison with the match-counting
  method and MCMC convergence diagnostics.
- **Estimate PCR-corrected drug efficacy** per site from the classification
  results and the number of patients enrolled.
- **Download** the posterior probabilities, the match-counting comparison and
  the convergence summary as CSV files.

To keep the shared server responsive, an analysis in the dashboard is limited
to **200 patients with a recurrence**, **1,000 additional Day 0 samples**,
**3 MCMC chains** and **5,000 iterations per chain**. Larger studies can be
split by site and analysed separately, or analysed with the MalReBay R package
below, which has no such limits.

## The MalReBay R package

This dashboard is built on MalReBay, an open-source R package. Using the
package directly provides full control over the analysis and supports
reproducible, documented workflows.

- [Package website](https://swisstph.github.io/MalReBay): overview,
  installation instructions and the function reference.
- [Vignette](https://swisstph.github.io/MalReBay/articles/MalReBay.html): a
  complete walkthrough of an analysis, from data import to the interpretation
  of results.
- [Analysis notebooks](https://github.com/SwissTPH/MalReBay/tree/main/inst/notebooks):
  focused notebooks for the analysis of microsatellite, MSP1/MSP2/GLURP and
  amplicon sequencing data.
- [Input data templates](https://github.com/SwissTPH/MalReBay/tree/main/inst/templates):
  Excel and CSV templates for each data type.
- [Source code](https://github.com/SwissTPH/MalReBay): the package repository
  on GitHub, for issue reports and contributions.

## Background

Therapeutic efficacy studies follow patients for several weeks after treatment
to assess the efficacy of antimalarial medicines. When parasites reappear
during follow-up, only a recrudescence is counted as a treatment failure. The
distinction is made by comparing the alleles detected in the parasite
population before and after treatment, and the resulting efficacy estimates
inform national treatment policies.

## Approach

MalReBay implements a Bayesian approach to this classification, which
complements the established allele match-counting method. For each patient, it
estimates the posterior probability that the recurrence is a recrudescence.
The model accounts for the local frequency of each allele, the similarity
between alleles, infections with multiple parasite clones and missing
genotyping data. It supports microsatellite, MSP1/MSP2/GLURP and amplicon
sequencing data.

This dashboard provides access to the method without running the package in R
or requiring technical expertise. The **How It Works** page of the dashboard
illustrates the underlying reasoning with an interactive example.

## Who we are

MalReBay was developed at Swiss TPH in Basel, Switzerland.

- **Monica Golumbeanu**, Swiss TPH, Basel, Switzerland
- **Veronica Adhiambo Ochieng**, Swiss TPH, Basel, Switzerland; African
  Institute for Mathematical Sciences (AIMS), Kigali, Rwanda

<a href="https://www.swisstph.ch"><img src="www/logos/SwissTPH_Logo_Screen_Indigo_rgb.png" alt="Swiss Tropical and Public Health Institute" height="60"></a>

## Funding

The development of MalReBay is funded by the Gates Foundation.

<a href="https://www.gatesfoundation.org"><img src="www/logos/Gates_Foundation_Logo.svg" alt="Gates Foundation" height="40"></a>
