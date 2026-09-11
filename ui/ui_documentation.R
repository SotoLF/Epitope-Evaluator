# ui_documentation.R -- input formats and supported predictors

.ee_pred_card <- function(name, url, paper, version, blurb, cols,
                          file = NULL, sep = "\t", header_rows = 1L) {
  bslib::accordion_panel(
    name,
    HTML(sprintf("<p>%s</p>", blurb)),
    HTML(sprintf(
      "<ul class='ee-meta'>
         <li><b>Version:</b> %s</li>
         <li><b>Paper:</b> %s</li>
         <li><b>Site:</b> <a href='%s' target='_blank' rel='noopener'>%s</a></li>
         <li><b>Columns used:</b> %s</li>
       </ul>", version, paper, url, url, cols)),
    if (!is.null(file)) ee_file_preview(file, sep = sep, header_rows = header_rows)
  )
}

documentation_ui <- function() {
  div(
    class = "ee-page",

    bslib::card(
      bslib::card_header("Input files"),
      bslib::card_body(HTML(
        "<p>Two files are required.</p>
         <p><b>1. The prediction file</b> &mdash; the output of a supported predictor, exactly
         as it was downloaded. Do not open and re-save it in a spreadsheet: NetMHC's
         <code>.xls</code> files are tab-separated text, and Excel will rewrite the two-row
         header that identifies the alleles.</p>
         <p><b>2. The FASTA file</b> &mdash; the multi-FASTA submitted to the predictor. It is
         needed for two things the prediction files do not reliably contain: the full protein
         identifiers (predictors truncate them to a fixed width) and the true protein lengths.</p>
         <p>The predictor is detected automatically from the file's header; you can override
         the guess. MHC class is inferred from the allele names, and sets the default cutoffs
         (2&nbsp;%rank for class&nbsp;I, 10 for class&nbsp;II).</p>
         <p><b>Score type.</b> A <b>percentile rank</b> is the rank of a peptide's predicted
         binding against a background of random natural peptides &mdash; lower is better, and
         it is comparable across alleles. A <b>score</b> is the predictor's raw output, either
         a binding affinity in nM or an elution-ligand score. Both work throughout the app;
         only the sensible cutoffs differ.</p>"),
        ee_figure("ui_input_upload.png",
                  "The upload panel. The predictor is detected from the file header and pre-selected; the choice can still be overridden.",
                  "?page=Analyse&tool=Input+data"))
    ),

    bslib::card(
      bslib::card_header("What a loaded dataset looks like"),
      bslib::card_body(
        HTML("<p>Whatever the source format, everything is normalised into one table:
              <code>Peptide</code>, <code>Pos</code>, <code>End</code>, peptide length,
              protein, protein length, then one column per MHC allele holding that
              peptide's score. Every tool is a view over this one table, and the
              <b>Download parsed table</b> button writes it in full.</p>
              <p>The summary strip reports what was read and which MHC class was detected.
              Anything the parser had to work around appears there as a note.</p>"),
        ee_figure("ui_input_loaded.png",
                  "A loaded dataset: 14,303 peptide/protein rows over 17 proteins and 12 class I alleles.",
                  "?page=Run+example&tool=Input+data"))
    ),

    bslib::card(
      bslib::card_header("Supported predictors"),
      bslib::card_body(bslib::accordion(
        open = FALSE,
        .ee_pred_card(
          "NetMHCpan", "https://services.healthtech.dtu.dk/service.php?NetMHCpan-4.1",
          "Reynisson et al. (2020). NetMHCpan-4.1 and NetMHCIIpan-4.0: improved predictions of MHC antigen presentation by concurrent motif deconvolution and integration of MS MHC eluted ligand data. <i>Nucleic Acids Research</i> 48:W449.",
          "4.1",
          "Pan-specific MHC class I predictor. Export the tabular output with
           <code>-xls -xlsfile out.xls</code>. Note the <b>two header rows</b>: the first
           carries an allele name at the start of each column block, the second names the
           columns within it. That is how the app finds the blocks, so it works whether or
           not the binding-affinity (<code>-BA</code>) columns are present. Re-saving through
           a spreadsheet destroys this structure.",
          "<code>Pos</code>, <code>Peptide</code>, <code>ID</code>, then <code>EL_Rank</code> (Rank) or <code>EL-score</code> (Score) per allele",
          file = "data/example_NetMHCPAN.xls", header_rows = 2L),

        .ee_pred_card(
          "NetMHCIIpan", "https://services.healthtech.dtu.dk/service.php?NetMHCIIpan-4.0",
          "Reynisson et al. (2020), as above.", "4.0",
          "Pan-specific MHC class II predictor. Same two-row header, plus a
           <code>Target</code> metadata column and three columns per allele. Positions are
           1-based here and 0-based in the class I tools; both are detected automatically.",
          "<code>Pos</code>, <code>Peptide</code>, <code>ID</code>, then <code>Rank</code> or <code>Score</code> per allele",
          file = "data/example_NetMHCIIPAN.xls", header_rows = 2L),

        .ee_pred_card(
          "NetMHC", "https://services.healthtech.dtu.dk/service.php?NetMHC-4.0",
          "Andreatta &amp; Nielsen (2016). Gapped sequence alignment using artificial neural networks: application to the MHC class I system. <i>Bioinformatics</i> 32:511.",
          "4.0",
          "Allele-specific MHC class I predictor. Three columns per allele, and its Score
           column is a binding affinity in <b>nM</b> &mdash; lower is stronger, and the usual
           threshold is 500 nM rather than a percentile.",
          "<code>Pos</code>, <code>Peptide</code>, <code>ID</code>, then <code>Rank</code> or <code>nM</code> per allele",
          file = "data/example_NetMHC.xls", header_rows = 2L),

        .ee_pred_card(
          "MHCflurry", "https://github.com/openvax/mhcflurry",
          "O&#39;Donnell et al. (2020). MHCflurry 2.0: improved pan-allele prediction of MHC class I-presented peptides by incorporating antigen processing. <i>Cell Systems</i> 11:42.",
          "2.0",
          "Comma-separated, and <b>long</b> rather than wide: one row per peptide and sample,
           naming the winning allele in <code>best_allele</code>. The app pivots this into a
           peptide &times; allele matrix. Pairs the predictor never scored stay missing rather
           than being read as a score of zero.",
          "<code>sequence_name</code>, <code>pos</code>, <code>peptide</code>, <code>best_allele</code>, then <code>affinity_percentile</code> (Rank) or <code>affinity</code> (Score)",
          file = "data/example_MHCFlurry.txt", sep = ",", header_rows = 1L),

        .ee_pred_card(
          "IEDB consensus", "http://tools.iedb.org/mhci/",
          "Moutaftsi et al. (2006). A consensus epitope prediction approach identifies the breadth of murine T(CD8+)-cell responses to vaccinia virus. <i>Nature Biotechnology</i> 24:817.",
          "2.24",
          "Tab-separated and long. <code>seq_num</code> is a <b>1-based index into the
           submitted FASTA</b>, not a name, so the FASTA must be the one that was uploaded and
           in the same order. Missing values are written as <code>-</code> and read as missing.",
          "<code>allele</code>, <code>seq_num</code>, <code>start</code>, <code>peptide</code>, then <code>consensus_percentile_rank</code> (Rank) or <code>ann_ic50</code> (Score)",
          file = "data/example_IEDB_consensus.txt", header_rows = 1L),

        bslib::accordion_panel(
          "Other (generic table)",
          HTML("<p>Any tab- or comma-separated table with a peptide, a start position, a
                protein identifier, a protein length, and then one column per MHC allele.
                Column names are used when they are recognisable (<code>peptide</code>,
                <code>pos</code>/<code>start</code>, <code>id</code>, <code>length</code>);
                otherwise that order is assumed. Every remaining column is treated as an
                allele, and its name becomes the label used throughout the app.</p>
               <p>A FASTA is optional for this format &mdash; without one, protein lengths come
                from the length column.</p>
               <p>The example below is the class I dataset from the paper. Its allele names
                arrive in the dotted form <code>HLA.A01.01</code>, which earlier versions of
                this app produced; they are canonicalised back to <code>HLA-A01:01</code> on
                load so that MHC class is detected and the labels read correctly.</p>"),
          ee_file_preview("Biological_Applications_Data/SARS_Class_I.txt", header_rows = 1L)
        )
      ))
    ),

    bslib::card(
      bslib::card_header("Notes on interpretation"),
      bslib::card_body(HTML(
        "<ul>
           <li><b>Unique peptides vs rows.</b> A peptide occurring in several proteins appears
               once per protein in the input. The Distribution, Intersection, Promiscuity and
               Conservation tools count <i>unique peptides</i>; Density counts occurrences,
               because the same sequence in two proteins really is two epitopes to present.</li>
           <li><b>Intersection vs Union.</b> Intersection takes the worst rank across the
               selected alleles, so a peptide qualifies only if it binds them all &mdash; the
               epitopes a whole population could present. Union takes the best rank, so binding
               any one allele is enough &mdash; what a single heterozygous individual sees.</li>
           <li><b>Missing values.</b> MHCflurry and IEDB do not score every peptide against
               every allele. Those pairs are treated as missing, not as non-binders, and are
               excluded from the min/max rather than dragging the whole row to NA.</li>
           <li><b>Density</b> is epitopes per amino acid, using the true FASTA length.</li>
         </ul>"))
    )
  )
}
