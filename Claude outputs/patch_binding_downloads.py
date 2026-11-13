#!/usr/bin/env python3
"""Binding > Similar Profiles: separate bar-plot / density-plot downloads with figure-friendly sizes.
Layout: 2x2 grid (bar plots on top, density plots below, each 8 in wide) + pooled density (16 in wide) underneath.
Run from the CHARM repo folder:   python3 patch_binding_downloads.py
Edits app.R and helper_functions.R in place (backups saved as *.bak). Safe to run twice."""
import re, shutil, sys

def load(p): return open(p, encoding="utf-8").read()
def save(p, s):
    shutil.copy(p, p + ".bak"); open(p, "w", encoding="utf-8").write(s)

# ---------------- helper_functions.R ----------------
h = load("helper_functions.R")
if "bar_plot_only" in h:
    print("helper_functions.R: already patched")
else:
    old = '  tbl_cols <- c("Profile", "RBP", "Target", "Correlation", if (has_ref) c("Percentile", "SameTarget"))\n  if (has_ref) {'
    assert h.count(old) == 1, "anchor 1 not found in helper_functions.R"
    h = h.replace(old, '  tbl_cols <- c("Profile", "RBP", "Target", "Correlation", if (has_ref) c("Percentile", "SameTarget"))\n'
                       '  bar_plot_only <- heatmap_plot   # bar plot alone (separate download)\n'
                       '  dens_only     <- NULL           # density plot alone (separate download)\n'
                       '  if (has_ref) {')
    old = '    heatmap_plot <- ggpubr::ggarrange(heatmap_plot, dens, ncol = 1, heights = c(1, 0.9))\n  }\n  list(\n    correlation_table = cor_df[, tbl_cols],\n    heatmap           = heatmap_plot\n  )'
    assert h.count(old) == 1, "anchor 2 not found in helper_functions.R"
    h = h.replace(old, '    bar_plot_only <- heatmap_plot\n    dens_only     <- dens\n'
                       '    heatmap_plot <- ggpubr::ggarrange(heatmap_plot, dens, ncol = 1, heights = c(1, 0.9))\n  }\n'
                       '  list(\n    correlation_table = cor_df[, tbl_cols],\n    heatmap           = heatmap_plot,\n'
                       '    bar_plot          = bar_plot_only,\n    density_plot      = dens_only\n  )')
    save("helper_functions.R", h); print("helper_functions.R: patched")

# ---------------- app.R ----------------
a = load("app.R")
old_v1 = '"_", dir_local, "_density"),\n            width = 16, height = 4'
if old_v1 in a:      # first version of this script was run: only the density width changes
    a = a.replace(old_v1, '"_", dir_local, "_density"),\n            width = 8, height = 4')
    save("app.R", a); print("app.R: density width corrected (16 -> 8)"); sys.exit(0)
if 'dir_local, "_barplot"' in a:
    print("app.R: already patched"); sys.exit(0)

# 1. density downloads: full figure width (2 x 8 in) x 4 in
a = a.replace('paste0(prefix, "_density"), width = 8, height = 4)', 'paste0(prefix, "_density"), width = 16, height = 4)')

# 2. store separately
old = "            dyn_dl_store[[plot_id]] <- res$heatmap\n            res$heatmap"
assert a.count(old) == 1, "anchor store not found"
a = a.replace(old, "            corr_store_plots(plot_id, res)\n            res$heatmap")

# 3. handlers for the 6 bar+density plots
old = '''          dl_id <- paste0("dl_", plot_id)
          output[[dl_id]] <- make_dyn_dl_handler(
            plot_id,
            paste0("binding_similar_profile_", ds_key, "_", dir_local)
          )'''
assert a.count(old) == 1, "anchor handler not found"
a = a.replace(old, '''          # bar plot and density plot are downloaded as two independent images
          output[[paste0("dl_", plot_id, "_bar")]] <- make_dyn_dl_handler(
            paste0(plot_id, "_bar"),
            paste0("binding_similar_profile_", ds_key, "_", dir_local, "_barplot"),
            width = 8, height = 6
          )
          output[[paste0("dl_", plot_id, "_density")]] <- make_dyn_dl_handler(
            paste0(plot_id, "_density"),
            paste0("binding_similar_profile_", ds_key, "_", dir_local, "_density"),
            width = 8, height = 4
          )''')

# 4. pooled (inc + dec) density: full width too
old = 'plot_id, paste0("binding_similar_profile_", ds_key, "_together")\n          )'
assert a.count(old) == 1, "anchor pooled not found"
a = a.replace(old, 'plot_id, paste0("binding_similar_profile_", ds_key, "_together"),\n            width = 16, height = 4\n          )')

# 5. buttons
pat = re.compile(r'downloadButton\("dl_binding_sim_plot_(both|K562|HEPG2)_(inc|dec)", "[^"]*Download plot",\n\s*class = "btn btn-sm btn-default"\)')
assert len(pat.findall(a)) == 6, "expected 6 buttons, found %d" % len(pat.findall(a))
a = pat.sub(lambda m: ('tagList(downloadButton("dl_binding_sim_plot_%s_%s_bar", "\\u2193 Download bar plot", '
                       'class = "btn btn-sm btn-default", style = "margin-right:6px;"), '
                       'downloadButton("dl_binding_sim_plot_%s_%s_density", "\\u2193 Download density plot", '
                       'class = "btn btn-sm btn-default"))') % (m.group(1), m.group(2), m.group(1), m.group(2)), a)
save("app.R", a); print("app.R: patched")
