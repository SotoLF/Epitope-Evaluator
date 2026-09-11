# ui_about.R -- the About page
#
# Static content only, so unlike v1 there is no about_server(): text that never
# changes does not need a renderUI() round-trip on every session start.

about_ui <- function() {
  div(
    class = "ee-page",

    bslib::card(
      bslib::card_header("What Epitope-Evaluator is"),
      bslib::card_body(HTML(
        "<p class='lead'>Epitope-Evaluator turns the output of a T-cell epitope predictor into
         six interactive analyses, so you can go from a table of scores to the epitopes worth
         following up without writing any code.</p>
         <p>Predictors such as NetMHCpan report a binding score for every peptide against every
         allele you asked for. For a whole proteome and a realistic HLA panel that is easily a
         million numbers &mdash; enough to answer almost any question, and far too many to read.
         Epitope-Evaluator applies one cutoff consistently across six views of the same table:
         how many epitopes there are, which alleles share them, which proteins carry them,
         where along a protein they sit, which bind the most alleles, and which survive across
         strains.</p>"),
        div(class = "ee-hero",
            ee_figure("tool_viewer.png",
                      "Epitope Viewer: predicted epitopes placed along the SARS-CoV-2 membrane protein, coloured by how many MHC alleles bind each one.",
                      "?page=Run+example&tool=Viewer")))
    ),

    bslib::layout_columns(
      col_widths = c(6, 6),
      bslib::card(
        bslib::card_header("What you need"),
        bslib::card_body(HTML(
          "<ol>
             <li>A <b>prediction file</b> from NetMHC, NetMHCpan, NetMHCIIpan, MHCflurry or
                 IEDB consensus &mdash; or any table in the documented generic layout.</li>
             <li>The <b>FASTA file</b> that was submitted to the predictor. It supplies the
                 real protein identifiers and lengths, which the prediction files truncate
                 or omit.</li>
           </ol>
           <p>Choose whether the scores are <b>percentile ranks</b> or raw <b>binding
           affinity / elution scores</b>. MHC class is detected from the allele names, and
           the default cutoffs follow it: 2&nbsp;%rank for class&nbsp;I, 10 for class&nbsp;II.</p>
           <p>No data leaves the session. Uploaded files live in a temporary directory for
           the life of the session and are discarded when it ends.</p>"))
      ),
      bslib::card(
        bslib::card_header("The six tools"),
        bslib::card_body(HTML(
          "<dl class='ee-dl-list'>
             <dt>Distribution</dt><dd>How epitopes spread across the score range, per allele or across allele combinations.</dd>
             <dt>Intersection</dt><dd>Which epitopes are shared between allele combinations, and which are private to one.</dd>
             <dt>Density</dt><dd>Epitope count against protein length; which proteins are unusually epitope-rich.</dd>
             <dt>Viewer</dt><dd>Where the epitopes sit along a single protein.</dd>
             <dt>Promiscuity</dt><dd>Epitopes binding many alleles at once, split into strong and weak binders.</dd>
             <dt>Conservation</dt><dd>Which epitopes survive across proteins, strains or variants.</dd>
           </dl>
           <p class='ee-hint'>Each tool links to a live example in the walkthrough on the
           <b>Tutorial</b> tab.</p>"))
      )
    ),

    bslib::card(
      bslib::card_header("What the tools look like"),
      bslib::card_body(
        bslib::layout_columns(
          col_widths = c(6, 6),
          ee_figure("tool_distribution.png", "Distribution \u2014 epitopes across the score range.",
                    "?page=Run+example&tool=Distribution"),
          ee_figure("tool_intersection.png", "Intersection \u2014 epitopes shared between alleles.",
                    "?page=Run+example&tool=Intersection"),
          ee_figure("tool_density.png", "Density \u2014 epitope count against protein length.",
                    "?page=Run+example&tool=Density"),
          ee_figure("tool_promiscuity.png", "Promiscuity \u2014 the broadest-binding epitopes.",
                    "?page=Run+example&tool=Promiscuity")
        ))
    ),

    bslib::card(
      bslib::card_header("Citation and contact"),
      bslib::card_body(HTML(
        "<p>If Epitope-Evaluator contributes to your work, please cite:</p>
         <blockquote class='ee-cite'>Soto, L. F., Requena, D., &amp; Fuxman Bass, J. I. (2022).
         Epitope-Evaluator: An interactive web application to study predicted T-cell epitopes.
         <i>PLoS ONE</i>, 17(8), e0273577.
         <a href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9417011/' target='_blank' rel='noopener'>PMC9417011</a>
         </blockquote>
         <p>Juan Fuxman Bass &mdash; fuxman@bu.edu<br>
            Luis F. Soto &mdash; lufesu98@gmail.com</p>
         <p>Source code and issue tracker:
            <a href='https://github.com/SotoLF/Epitope-Evaluator' target='_blank' rel='noopener'>github.com/SotoLF/Epitope-Evaluator</a>.
            Released under the MIT licence.</p>")),
      bslib::card_footer(class = "ee-version",
                         sprintf("Version %s · R %s", EE$version, getRversion()))
    )
  )
}
