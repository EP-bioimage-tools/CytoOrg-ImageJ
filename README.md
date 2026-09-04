# CytoOrg

**CytoOrg** is a Fiji/ Groovy script for quantifying cytoskeletal organization in 2D fluorescence microscopy images. It analyzes one cytoskeletal signal at a time and supports both isolated single-cell images and fields containing multiple cells described by binary masks.

## What it measures

The pipeline enhances filamentous structures using a Tubeness filter, extracts their **skeleton**, and quantifies:

- cell geometry (area, roundness, solidity);
- filament thickness, length, density, tortuosity and branching;
- global and local orientation statistics, including angular distributions and n. of direction changes per µm.

## Requirements

- Current **Fiji** installation with Fiji updates applied;
- **AnalyzeSkeleton**, **Local Thickness**, and **Tubeness** (pre-installed in a standard Fiji installation);
- **OrientationJ**, installed through the Fiji updater (BIG-EPFL update site) or as a Fiji plugin.

## Installation and run

1. Download [`CytoOrg_v1.0.groovy`](CytoOrg_v1.0.groovy).
2. Open it in Fiji's **Script Editor**.
3. Select **Groovy** as the language and click **Run**.
4. Select the input folder, choose the analysis options, and start the analysis.

Results are written to a `Results` folder inside the selected input folder by default.

## Input images

- Supported formats: `.tif`, `.tiff`, `.png`, `.jpg`, and `.jpeg` (case-insensitive).
- 2D single planes or maximum-intensity projections; one cytoskeletal channel per image.
- The complete cell should be visible and the pixel size should be known.
- Input files must be directly inside the selected folder; subfolders are not scanned.
- The input folder must be writable.
- Files beginning with `segcell_` are reserved for multi-cell masks and are not analyzed as primary images.

### Multi-cell masks (optional)

Enable **Use segcell masks for multiple cells** and place an 8-bit binary mask next to each image. The mask must have identical dimensions and be named:

```text
segcell_<complete original filename>
```

## Output

The main summary is `Results.csv`. Depending on the selected options, the script can also save:

- thickness, skeleton, branching, orientation, and angular-distribution .tif maps;
- optional enhanced maps with thicker lines for visualization;
- individual- and grouped-filament CSV tables containing the raw data of AnalyzeSkeleton;

The repository includes colour bar/key images for interpreting the thickness, branching, and orientation maps. 

For complete instructions, input requirements, parameter descriptions, and troubleshooting, see the [CytoOrg Guide v1.0](<CytoOrg Guide v1.0.pdf>).


## Citation, license, and contact

CytoOrg is released under the **MIT License** (© 2026 Edoardo Perina).

For questions or bug reports: **e.perina@studenti.unibs.it**

Repository: <https://github.com/EP-bioimage-tools/CytoOrg-ImageJ>

Images shown in the guide are from the IDR0050 experimentA dataset, licensed under **CC BY 4.0** (DOI: <https://doi.org/10.17867/10000119>).
