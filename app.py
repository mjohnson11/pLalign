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
from Bio import SeqIO
import mappy as mp
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit.components.v1 as components
from string import Template
import json

annotation_length_thresh = 50
alignment_len_thresh = 10
max_alignments = 1000
mapq_thresh = 0

THRESH = 1000

def get_genbank_features(ref_data):
    ref_features = []
    for f in ref_data.features:
        start = f.location.start
        end = f.location.end
        strand = f.location.strand
        label = f.qualifiers.get('label', [''])[0] if 'label' in f.qualifiers else ''
        note = f.qualifiers.get('note', [''])[0] if 'note' in f.qualifiers else ''
        if (end-start) > annotation_length_thresh:
            ref_features.append({'start': start, 'end': end, 'strand': strand, 'label': label, 'note': label + ' ' + note})
    return ref_features

def do_alignment(aligner, read_data, stop_point=None):
    i = 0
    all_hits = []
    read_lens = []
    for rec in read_data: #name, seq, qual in mp.fastx_read(io.StringIO(text)):
        if type(rec) == str:
            seq = rec
        else:
            seq = str(rec.seq)
        read_lens.append(len(seq))
        i += 1
        if i > THRESH:
            break
        # hits sorted by length
        hits = sorted([h for h in aligner.map(seq)], key=lambda hit: -1*(hit.r_en-hit.r_st))
        if len(hits) > 0:
            hit = hits[0]
            hit_len = hit.r_en-hit.r_st
            if hit_len > alignment_len_thresh and hit.mapq >= mapq_thresh:
                all_hits.append({
                        'cigar': hit.cigar, 
                        'ref_start': hit.r_st, 
                        'ref_end': hit.r_en, 
                        'strand': hit.strand, 
                        'query_start': hit.q_st, 
                        'query_end': hit.q_en, 
                        'query_seq': seq,
                        'mapq': hit.mapq
                    })

        if (stop_point and stop_point < i) or len(all_hits)>max_alignments:
            break
    st.write(len(all_hits))
    return all_hits, read_lens
        


# Title and description
st.title("Plasmid Alignment App")
st.write("Align reads to a reference plasmid and visualize the results.")

# Radio button for reference input
ref_type = st.radio("Reference Type", ("File Upload", "Pasted Sequence", "Test data"))
ref_seq = None
if ref_type == "Test data":
    ref_data = SeqIO.read('./test_data/pmj09.gb', "genbank")
    ref_seq = ref_data.seq
    ref_features = get_genbank_features(ref_data)
elif ref_type == "File Upload":
    # Upload reference file
    ref_file = st.file_uploader("Upload Reference (GenBank or FASTA)", type=["gb", "gbk", "fasta"])
    if ref_file:
        # Parse reference sequence
        is_genbank = '.gb' in ref_file.name
        text_io = io.TextIOWrapper(ref_file, encoding="UTF-8")
        text = text_io.read()
        st.success("File uploaded.")
        ref_data = SeqIO.read(io.StringIO(text), "genbank" if is_genbank else "fasta")
        ref_seq = ref_data.seq
        ref_features = []
        if is_genbank:
            ref_features = get_genbank_features(ref_data)
            
        st.write(len(ref_features))
else:
    # Text area for pasted sequence
    ref_seq = st.text_area("Paste Reference Sequence")
    ref_features = []
    
# Validate reference sequence
if not ref_seq:
    st.warning("Please provide a valid reference sequence.")
    st.stop()
    
all_hits = None
# Upload reads file
if ref_type == 'Test data':
    aligner = mp.Aligner(seq=str(ref_seq)+str(ref_seq)) # DOUBLING REF SEQ because it's circular
    read_data = SeqIO.parse('./test_data/FAY32610_pass_barcode13_92612969_9c09aa16_0.fastq', "fastq")
    all_hits, read_lens = do_alignment(aligner, read_data)
else:
    aligner = mp.Aligner(seq=str(ref_seq)+str(ref_seq)) # DOUBLING REF SEQ because it's circular
    reads_type = st.radio("Reads Type", ("File Upload", "Pasted Sequence"))
    if reads_type == "File Upload":
        reads_file = st.file_uploader("Upload Reads (FASTQ or FASTA)", type=["fastq", "fasta"])
        if reads_file:
            # Minimap2 alignment
            text_io = io.TextIOWrapper(reads_file, encoding="UTF-8")
            text = text_io.read()
            st.success("File uploaded.")
            read_data = SeqIO.parse(io.StringIO(text), "fastq" if reads_file.name.endswith(".fastq") else "fasta")
            all_hits, read_lens = do_alignment(aligner, read_data)
    else:
        read_seq = st.text_area("Paste Reference Sequence")
        if read_seq:
            all_hits, read_lens = do_alignment(aligner, [read_seq])

if all_hits:
    #process_alignments(all_hits, str(ref_seq))
    
    # Placeholder for visualization
    st.write("Alignment Results:")
    
    # 2D structure for odd cases with more than one alignment per read
    alignment_lens = [h['ref_end']-h['ref_start'] for h in all_hits]

    # Create traces for each histogram
    trace1 = go.Histogram(x=read_lens, name='Read Lengths', opacity=0.6)
    trace2 = go.Histogram(x=alignment_lens, name='Alignment Lengths', opacity=0.6)

    # Create the figure with stacked histograms
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True)
    fig.add_trace(trace1,row=1, col=1)
    fig.add_trace(trace2,row=2, col=1)
  
    # Display the figure in Streamlit
    st.plotly_chart(fig)
    
    full_data = json.dumps({
        'ref_seq': str(ref_seq),
        'ref_features': ref_features,
        'alignments': all_hits
    })

    # ... (JavaScript and d3 code for visualization)
    html_template = Template(
        """
        <div id="pLalignViz" style="width:100%"></div>
        <style>
        $pLalignCSS
        </style>
        <script src="https://d3js.org/d3.v6.min.js"></script>
        <script>
        const pLalignData = $pLalignData
        $pLalignScript
        console.log(pLalignData)
        </script> 
        """
    )
    js_text = ''
    # for offline dev
    #with open('./d3.v6.min.js', 'r') as infile:
    #    js_text += infile.read()
    with open('./pLalign.js', 'r') as infile:
        js_text += infile.read()
        
    with open('./pLalign.css', 'r') as infile:
        css_text = infile.read()
    components.html(html_template.substitute(pLalignScript=js_text, pLalignData=full_data, pLalignCSS=css_text), height=1400, width=1000)
else:
    st.warning("Please provide a reads file.")
    st.stop()
    
    
    
jnk  = """
## OLD function from starting to refactor for python processing - broken somehow, caused app to break

def refpos(pos, refLen):
    while pos > refLen-1:
      pos = pos - refLen
    if pos < 0:
        return refLen + pos
    else:
        return pos

def process_alignments(als, refSeq):
    refLen = len(refSeq)
    c = 0
    for a in als:
        c += 1
        print(c)
        query_position = a['query_start']
        ref_position = a['ref_start']
        a['mismatches'] = []
        a['deletions'] = []
        a['insertions'] = []
        a['alignment_by_ref_pos'] = [[]]*refLen
        if query_position > 0:
            for i in range(query_position):
                a['alignment_by_ref_pos'][refpos(ref_position-query_position+i, refLen)].append(['pre-alignment', a['query_seq'][i]])
        for cig in a['cigar']:
            bp = cig[0]
            letter = "MIDNSHP=XB"[cig[1]]
            if letter == 'S':
                for i in range(bp):
                    a['alignment_by_ref_pos'][refpos(ref_position+i, refLen)].append(['pre-alignment', '+', a['query_seq'][query_position+i]])
                query_position += bp
            elif letter == 'M':
                for i in range(bp):
                    if a['query_seq'][query_position+i] != (refSeq+refSeq)[ref_position+i]:
                        a['mismatches'].append(ref_position+i)
                        a['alignment_by_ref_pos'][refpos(ref_position+i, refLen)].append(['mismatch', a['query_seq'][query_position+i]])
                    else:
                        a['alignment_by_ref_pos'][refpos(ref_position+i, refLen)].append(['match', a['query_seq'][query_position+i]])
                query_position += bp
                ref_position += bp
            elif letter == 'D':
                for i in range(bp):
                    a['alignment_by_ref_pos'][refpos(ref_position+i, refLen)].append(['deletion', '-'])
                ref_position += bp
                a['deletions'].append([ref_position-bp, ref_position])
            elif letter == 'I':
                a['alignment_by_ref_pos'][refpos(ref_position, refLen)].append(['insertion', '', a['query_seq'][query_position:query_position+bp]])
                query_position += bp
                a['insertions'].append([ref_position, bp])

    if query_position < len(a['query_seq'])-1:
        print('Q', len(a['query_seq']), query_position)
        for i in range(len(a['query_seq'])-1-query_position):
            a['alignment_by_ref_pos'][refpos(ref_position+i, refLen)].append(['post-alignment', a['query_seq'][query_position+i]])

""";