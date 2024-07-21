### TODO

# refactor so alignment strings etc. are made in python
# ref features in inset
# insertion rectangle showing extent in inset
# controls for zoom and drag in inset
# options for thresholds at top
# options for bp start - end for barcode extraction and clustering after viz.
# bigger shot - linear viz option


import streamlit as st
import io
import numpy as np
from Bio import SeqIO

# Title and description
st.title("Plasmid Alignment App")
st.write("Align reads to a reference plasmid and visualize the results.")
st.write(str(np.exp(4)))