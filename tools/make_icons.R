# tools/make_icons.R -- generate the app mark and its favicon fallbacks
#
#   Rscript tools/make_icons.R
#
# The mark echoes the Epitope Viewer, which is the app's most distinctive
# figure: a protein backbone with epitope blocks beneath it, coloured on the
# same yellow-to-red "number of alleles bound" ramp. Drawn from the palette in
# utils/constants.R so the branding and the plots cannot drift apart.
#
# SVG is written by hand (Chrome, Firefox and Safari all take an SVG favicon);
# the PNGs are drawn with grDevices for older browsers and iOS, which avoids
# depending on rsvg/ImageMagick being installed.

source("utils/constants.R")

NAVY    <- "#0B3C5D"
PROTEIN <- EE$col$protein          # #9CC3E5
EP      <- c("#FFCE5C", "#F08C3C", "#B2182B")

# Geometry on a 32x32 grid: x, y, width, height, fill.
# One backbone bar, then two rows of epitope blocks offset like a real layout.
BARS <- list(
  list(x =  5, y =  7.0, w = 22, h = 4.6, fill = PROTEIN),
  list(x =  5, y = 14.2, w = 10, h = 4.6, fill = EP[1]),
  list(x = 17, y = 14.2, w = 10, h = 4.6, fill = EP[2]),
  list(x =  9, y = 21.4, w =  9, h = 4.6, fill = EP[2]),
  list(x = 20, y = 21.4, w =  7, h = 4.6, fill = EP[3])
)

# ---------------------------------------------------------------------------
# SVG
# ---------------------------------------------------------------------------

svg_rect <- function(b, r = 2.3) {
  sprintf('  <rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" rx="%.1f" fill="%s"/>',
          b$x, b$y, b$w, b$h, r, b$fill)
}

write_svg <- function(path, bg = TRUE) {
  body <- c(
    # Explicit width/height as well as the viewBox: without an intrinsic size,
    # Firefox reports naturalWidth 0 and some contexts fail to size the image.
    '<svg xmlns="http://www.w3.org/2000/svg" width="32" height="32"',
    '     viewBox="0 0 32 32" role="img" aria-label="Epitope-Evaluator">',
    if (bg) sprintf('  <rect width="32" height="32" rx="7" fill="%s"/>', NAVY),
    vapply(BARS, svg_rect, character(1)),
    '</svg>'
  )
  writeLines(body, path)
  cat(sprintf("  %-26s %d bytes\n", path, file.size(path)))
}

# ---------------------------------------------------------------------------
# PNG fallbacks
# ---------------------------------------------------------------------------

rounded <- function(x, y, w, h, r, col) {
  # A rounded rectangle from two overlapping rects plus four corner circles.
  rect(x + r, y, x + w - r, y + h, col = col, border = NA)
  rect(x, y + r, x + w, y + h - r, col = col, border = NA)
  th <- seq(0, 2 * pi, length.out = 60)
  for (cx in c(x + r, x + w - r)) for (cy in c(y + r, y + h - r)) {
    polygon(cx + r * cos(th), cy + r * sin(th), col = col, border = NA)
  }
}

write_png <- function(path, px) {
  grDevices::png(path, width = px, height = px, bg = "transparent")
  op <- par(mar = rep(0, 4), xaxs = "i", yaxs = "i")
  # y is flipped so the SVG coordinates can be reused unchanged.
  plot.new(); plot.window(xlim = c(0, 32), ylim = c(32, 0))
  rounded(0, 0, 32, 32, 7, NAVY)
  for (b in BARS) rounded(b$x, b$y, b$w, b$h, 2.3, b$fill)
  par(op); grDevices::dev.off()
  cat(sprintf("  %-26s %dx%d, %d bytes\n", path, px, px, file.size(path)))
}

# ---------------------------------------------------------------------------

dir.create("www", showWarnings = FALSE)
write_svg("www/favicon.svg")
write_png("www/favicon-32.png", 32)
write_png("www/favicon-180.png", 180)
cat("\nDone. Referenced from the <head> block in app.R.\n")
