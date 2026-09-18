# Tutorial_Network_Analysis
+ `Cytoscape` (v3.10.4) was used for visualization.

+ `Experiment 1` in `Gut-Lung Axis` project
  (1) Indices from mouse *in vivo* experiment
  (2) Microbial read frequencies from mouse fecal shotgun metagenomics sequencing results (species-level)
      (results from 16S rRNA amplicon sequencing could be used instead)
  (3) Intensities of metabolic features from untargeted Orbitrap LC-MS/MS analysis

  above data were used for network analysis.


## Step 0 : Change raw data to Cytoscape format
### Co-occurrence
+ Use `05. correlation+cytoscape format.ipynb` script to make a Cytoscape format.

### Kernel causality
+ Use `Kernel Causality_GutLungAxis.R` script to make a Cytoscape format.


## Step 1 : Input
### Co-occurrence
1) `File` -> `Import` -> `Network from File...` -> `cytoscape_edges.csv`
2) `File` -> `Import` -> `Table from File...` -> `cytoscape_nodes.csv`

### Kernel causality
1) `File` -> `Import` -> `Network from File...` -> `edges_kernel.tsv`
2) `File` -> `Import` -> `Table from File...` -> `nodes_kernel.tsv`


## Step 2 : Visualization
### Co-occurrence

### Kernel causality
+ Edges
Label Font Size : 12
Line Type : Column-Polarity / Mapping Type-discrete
            pos-solid / neg-dashed
Stroke Color : #C9C9C9 (Edge color to arrows)
Target Arrow Shape : Column-Interaction / Mapping Type-discrete
                     causes-triangle
Width : Column-EdgeWidth / Mapping Type-passthrough

+ Nodes
Border Color : #C9C9C9
Border Width : 3
Fill Color : #FFB000
Height : 13 (Lock node width and height)
Label : Column-name / Mapping Type-passthrough
Label Font Size : 15
Label Position : bottom-right
Opacity : 80%

+ After visualization is done, save it to `png` format file by
  1) Click `Fit Content` button (its magnifier icon below `Help` button).
  2) `File` -> `Export` -> `Network to Image...` ->
  3) + Export File Format: `PNG (*.png)`
     + Save Image as: (to your output file path)
     + Zoom (%): `500%`


## Step 3 : Add legends
### Co-occurrence
+ This is the step to make legends for `Nodes` and `Edges`.
  + You can just modify `co-occurrence_legend_session_GutLungAxis` file to your own design in `Step 2`.
  + For `Edges`, you can use `[ppt양식] Cytoscape_coefficient+EdgeWidth_legends.pptx`.

