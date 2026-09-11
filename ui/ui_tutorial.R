# ui_tutorial.R -- walkthrough of the app
#
# Every screenshot here is a capture of THIS build, produced by
# tools/screenshots.py against the running app. Re-run that script after any
# layout change so the tutorial never shows a UI that no longer exists.
#
# v1 hard-coded eight <video> tags pointing at www/videos/*.mp4, which are
# gitignored and absent from the repository, so a fresh clone showed eight
# broken players. Videos are still used if present, and simply skipped if not.

.ee_step <- function(n, title, body, figures = list(), video = NULL) {
  has_video <- !is.null(video) && file.exists(file.path("www", "videos", video))
  bslib::accordion_panel(
    sprintf("%d. %s", n, title),
    HTML(body),
    lapply(figures, function(f) ee_figure(f$file, f$caption, f$href)),
    if (has_video) {
      tags$video(src = file.path("videos", video), type = "video/mp4",
                 controls = NA, preload = "none", class = "ee-video")
    }
  )
}

.lnk <- function(tool, page = "Run example") {
  sprintf("?page=%s&tool=%s", gsub(" ", "+", page), gsub(" ", "+", tool))
}

tutorial_ui <- function() {
  div(
    class = "ee-page",

    bslib::card(
      bslib::card_header("Getting started"),
      bslib::card_body(
        HTML(
          "<p>The quickest way in is the <b>Run example</b> tab: it loads a prediction for the
           SARS-CoV-2 proteome and every tool is immediately populated, with nothing to upload.
           Switching the example between class&nbsp;I and class&nbsp;II predictors is the fastest
           way to see the cutoffs follow the MHC class.</p>
           <p>Every tool works the same way. The sidebar is split in two:</p>
           <ul>
             <li><b>Above the Run button</b> &mdash; the parameters that define the analysis
                 (alleles, proteins, cutoffs, AND/OR). Changing them does nothing until you
                 press <b>Run analysis</b>, so a large recomputation never starts under you.</li>
             <li><b>Below, under DISPLAY</b> &mdash; presentation only (plot type, sort order,
                 log scale). These apply immediately.</li>
           </ul>
           <p>A tool tab is never blank: on arrival it already shows the analysis at the
           dataset's own defaults. Every figure saves as a PNG from the camera icon in its
           toolbar, and every table has a download button that writes the complete result,
           not just the rows on screen.</p>"),
        div(class = "ee-hero",
            ee_figure("tool_distribution.png",
                      "A tool page: parameters on the left (analysis above Run, display below), the summary strip, then the figures and tables.",
                      .lnk("Distribution")))
      )
    ),

    bslib::card(
      bslib::card_header("Walkthrough"),
      bslib::card_body(bslib::accordion(
        open = FALSE,

        .ee_step(1, "Load your data",
          "<p>On the <b>Analyse</b> tab, upload the predictor output and the FASTA it was run
           on. The format is detected from the file header and reported above the dropdown;
           override it if the guess is wrong. Choose whether the scores are percentile ranks
           or raw affinities, then press <b>Load data</b>.</p>
           <p>The summary strip reports how many peptides, proteins and alleles were read, and
           which MHC class was detected. Anything unusual about the input &mdash; protein IDs
           resolved by prefix matching, alleles dropped for having no numeric values, a FASTA
           whose lengths disagree with the prediction &mdash; appears there as a note rather
           than being silently absorbed.</p>",
          list(list(file = "ui_input_upload.png",
                    caption = "The upload panel. The predictor is auto-detected from the file header.",
                    href = .lnk("Input data", "Analyse")),
               list(file = "ui_input_loaded.png",
                    caption = "After loading: the summary strip, plus the parsed peptide x allele table you can download.",
                    href = .lnk("Input data"))),
          "UploadFiles.mp4"),

        .ee_step(2, "Distribution: pick a cutoff you can defend",
          "<p>Start here. The histogram shows where your peptides actually fall across the
           score range, and the cumulative view answers \"how many epitopes do I get if I move
           the threshold?\" directly. The heatmap underneath repeats the count for each allele
           separately across a ladder of cutoffs &mdash; the fastest way to see whether one
           allele is carrying your whole epitope set.</p>",
          list(list(file = "tool_distribution.png",
                    caption = "Three class I alleles under the Union condition: 442 epitopes at 2 %rank, and how that count grows with the cutoff.",
                    href = .lnk("Distribution"))),
          "EpitopeDistribution.mp4"),

        .ee_step(3, "Intersection: who shares what",
          "<p>Select the alleles of interest and read the UpSet plot: the bars on top are
           intersection sizes, the dot matrix says which alleles each bar refers to, and the
           grey bars on the left are each allele's total. A single dot is an epitope private
           to that allele; connected dots are epitopes shared by exactly those alleles. The
           large, highly-connected bars are your promiscuous candidates.</p>
           <p>For two to four alleles the Venn view is often easier to read; past that, UpSet
           is the only one that stays legible. Switching between them is a display option, so
           it applies without re-running.</p>",
          list(list(file = "tool_intersection.png",
                    caption = "Four alleles: 1,431 epitopes in total, 158 shared by two or more, in 12 distinct combinations.",
                    href = .lnk("Intersection")),
               list(file = "tool_intersection_venn.png",
                    caption = "The same question as a Venn diagram, drawn natively so no system graphics libraries are needed.",
                    href = .lnk("Intersection"))),
          "EpitopeIntersection.mp4"),

        .ee_step(4, "Density: find the epitope-rich proteins",
          "<p>Long proteins carry more epitopes simply because they contain more peptides, so
           look for points sitting well above the trend. Clicking rows in the table highlights
           them in the plot. The grid underneath breaks the same counts down by allele, which
           shows whether a protein is broadly immunogenic or presented by only a few HLA
           types.</p>",
          list(list(file = "tool_density.png",
                    caption = "17 proteins, 4,140 epitopes. ORF3C is the densest per amino acid, despite R1AB carrying the most in absolute terms.",
                    href = .lnk("Density"))),
          "EpitopeDensity.mp4"),

        .ee_step(5, "Viewer: where the epitopes are",
          "<p>Pick one protein and see its epitopes drawn to scale, stacked into rows so
           overlaps stay readable and coloured by how many alleles bind each one. Dense red
           clusters are the hotspots. Drag to zoom into a region; hover for the sequence,
           coordinates and allele list.</p>",
          list(list(file = "tool_viewer.png",
                    caption = "The membrane protein (222 aa): 76 epitopes covering 96% of the sequence, in 9 stacked rows.",
                    href = .lnk("Viewer"))),
          "EpitopeLocation.mp4"),

        .ee_step(6, "Promiscuity: the broadest binders",
          "<p>Set the strong and weak cutoffs and a minimum number of alleles. The heatmap
           lists epitopes sorted by how many alleles they bind, red for strong binders and
           orange for weak ones. The rows at the top are the candidates that would cover the
           most genetically diverse population.</p>",
          list(list(file = "tool_promiscuity.png",
                    caption = "Epitopes binding many of the twelve alleles at once, strong binders in red.",
                    href = .lnk("Promiscuity"))),
          "EpitopePromiscuity.mp4"),

        .ee_step(7, "Conservation: what survives across variants",
          "<p>Run the prediction on a multi-FASTA holding the same protein from several
           strains, then compare those proteins here. Epitopes in the all-proteins
           intersection are conserved targets; epitopes private to one variant mark positions
           where escape has already occurred.</p>
           <p>The example below compares four different SARS-CoV-2 proteins rather than
           variants of one, which is the same machinery applied to a different question:
           which epitopes recur across the proteome.</p>",
          list(list(file = "tool_conservation.png",
                    caption = "Four proteins compared; the overlapping regions are epitopes shared between them.",
                    href = .lnk("Conservation"))),
          "EpitopeConservation.mp4")
      ))
    ),

    bslib::card(
      bslib::card_header("Practical tips"),
      bslib::card_body(HTML(
        "<ul>
           <li><b>The tools are independent.</b> Each has its own Run button, so changing an
               allele selection in one does not recompute the others.</li>
           <li><b>Deep links work.</b> Any view can be shared as a URL, e.g.
               <code>?page=Run+example&amp;tool=Viewer</code>. The screenshots above link to
               the live view they illustrate.</li>
           <li><b>Downloads are complete.</b> Large tables are paged on screen, long peptide
               lists are shortened to keep rows readable, and the figures cap how much they
               draw &mdash; but every download button writes the full result.</li>
           <li><b>Figures export at print resolution.</b> The camera icon saves a
               1600&times;1000 PNG at 2&times; scale, not a screenshot of the panel.</li>
           <li><b>If a panel shows a red message</b>, that is the actual error from the
               analysis &mdash; usually a cutoff that excludes everything, or two cutoffs the
               wrong way round. The rest of the session keeps working.</li>
         </ul>"))
    )
  )
}
