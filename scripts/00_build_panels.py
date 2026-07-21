from voi_remsen import data, viz

def main():
      weekly = data.build_calhabmap_weekly()
      cdphw = data.build_cdph_weekly()
      cell = data.build_cell_panel()
      tox = data.build_toxicity_panel()
      joint = data.build_joint_panel()

      print(f"calhabmap_weekly: {len(weekly):5d} site-weeks, {weekly.location_code.nunique()} sites")
      print(f"cdph_weekly:      {len(cdphw):5d} site-weeks")
      print(f"cell_panel:       {len(cell):5d} weeks, {cell.location_code.nunique()} sites "
          f"(biomass; median {cell.cells.median():.0f} cells/L)")
      print(f"toxicity_panel:   {len(tox):5d} weeks; category counts "
          f"{tox.da_cat.value_counts().sort_index().to_dict()}; closures {int((tox.da>=20).sum())}")
      print(f"joint_panel:      {len(joint):5d} site-weeks; "
          f"cells-only {int((joint.cells.notna() & joint.da.isna()).sum())}, "
          f"DA-only {int((joint.da.notna() & joint.cells.isna()).sum())}, "
          f"both {int((joint.cells.notna() & joint.da.notna()).sum())}")

      viz.plot_coverage()
if __name__ == "__main__":
      main()
