# Paper Figure Style Guidelines

These notes collect figure-design lessons from the DStressR manuscript work and
are intended as reusable guidance for future Codex paper projects.

## Panel Labels

- Use small, bold, lowercase panel labels (`a`, `b`, `c`, ...) for manuscript
  figures unless the target journal requests another style.
- Place panel labels in the outer layout margin of the composite figure, near
  the top-left of each subpanel.
- Do not draw panel labels inside the data region.
- Do not write explicit labels such as `Figure 3a` inside the graphic. The
  manuscript caption and cross-references carry the figure number.
- For multi-panel figures composed from separate plots, add labels at the final
  composition stage rather than inside each ggplot panel.

## Captions and Text

- Keep titles and explanatory prose out of the plotting area. Use the caption
  to explain what each panel shows.
- Refer to panels as "panel a", "panels b--d", etc. in the caption and main
  text.
- The caption should define the displayed variables, estimands, and thresholds,
  especially when axes are method-specific.
- Avoid expansive top captions inside the figure. Facet strips or short column
  labels are acceptable when they identify methods or conditions.

## Scientific Plotting Defaults

- Start from the scientific question, not from the plot type. Each panel should
  show a distinct aspect of the analysis.
- Keep data regions clean. Avoid decorative labels, large legends, and
  annotation boxes that compete with the data.
- Use consistent palettes across figures when the same biological classes recur,
  for example promoter-specific colors.
- Use separate, clearly distinguishable colors for method-set diagrams such as
  Venn diagrams; avoid reusing the promoter palette for method categories.
- Check grayscale distinguishability when color carries method or group
  identity.
- Use point shapes as well as color when two classes must remain identifiable
  in print.

## Multi-Panel Layouts

- Match the layout to the reading order of the result: diagnostics first,
  primary effect plots next, then set-overlap or disagreement summaries.
- Align related panels by axis scale when direct comparison is intended.
- Use shared scales when the numerical comparison is meaningful, and explicitly
  state when scales are method-specific.
- Reserve large panels for dense label placement, such as volcano plots with
  annotated hits.
- Put legends inside panels only when they do not occlude data and when doing so
  saves figure space.

## P-Value and Hit-Calling Diagnostics

- Do not color entire histogram bins as "significant" when a bin contains both
  called and uncalled tests.
- P-value histograms should usually show raw p-values neutrally; describe FDR
  adjustment and hit rules in the caption or text.
- If significant calls are highlighted in p-value diagnostics, mark individual
  observations or use a separate summary plot rather than recoloring histogram
  bins.
- For method-comparison scatter plots, show all tests in a muted background and
  highlight only scientifically relevant discordant sets.

## Volcano Plots

- Use the method's native effect scale only when it is clearly stated in the
  caption.
- Keep the p-value axis comparable across panels when possible.
- Annotate biologically important or overlapping hits with thin leader lines.
- Spend time on label placement. A publication volcano plot should use open
  white space and avoid collisions rather than simply labeling the top points.

## Implementation Notes

- In R/ggplot workflows, build panel labels after arranging grobs, for example
  by overlaying `grid::textGrob()` labels on the final `gridExtra` or `patchwork`
  composition.
- Keep reusable helper functions close to the figure script until the style is
  stable; then promote them to shared plotting utilities.
- Save both PNG and PDF outputs for manuscript figures.
- Inspect rendered figures visually before trusting the manuscript render.
