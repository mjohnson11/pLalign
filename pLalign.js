
// PARAMETERS
const alignment_num = 30;
const background_color = '#333';
const base_color = 'white';
const colors = {
  'al': '#777',
  'X': '#F00',
  'I': '#0FF',
  'D': '#E4691D',
  'outside': '#333',
  'overlap': '#BBB',
  'button_text': '#FFF',
  'A': '#4665C9',
  'T': '#C2AA23',
  'G': '#5FB152',
  'C': '#B1362D',
  'N': '#777',
  'pre-alignment': 'green',
  'post-alignment': 'red'
}
colors['mismatch'] = colors['X'];
colors['match'] = colors['al'];
//colors['pre-alignment'] = colors['outside'];
//colors['post-alignment'] = colors['outside'];
colors['insertion'] = colors['I']; 
colors['deletion'] = colors['D'];

// UTILITY FUNCTIONS
function rc(s) {
  // reverse complement
  const complement = {'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 
                      'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', '-': '-'};
  let new_seq = [];
  for (let i=s.length-1; i>-1; i--) {
    new_seq.push(complement[s[i]] || s[i]);
  }
  return new_seq.join('').toUpperCase();
}

// https://gist.github.com/tophtucker/62f93a4658387bb61e4510c37e2e97cf
function measureText(string, fontSize = 10) {
  const widths = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.2796875,0.2765625,0.3546875,0.5546875,0.5546875,0.8890625,0.665625,0.190625,0.3328125,0.3328125,0.3890625,0.5828125,0.2765625,0.3328125,0.2765625,0.3015625,0.5546875,0.5546875,0.5546875,0.5546875,0.5546875,0.5546875,0.5546875,0.5546875,0.5546875,0.5546875,0.2765625,0.2765625,0.584375,0.5828125,0.584375,0.5546875,1.0140625,0.665625,0.665625,0.721875,0.721875,0.665625,0.609375,0.7765625,0.721875,0.2765625,0.5,0.665625,0.5546875,0.8328125,0.721875,0.7765625,0.665625,0.7765625,0.721875,0.665625,0.609375,0.721875,0.665625,0.94375,0.665625,0.665625,0.609375,0.2765625,0.3546875,0.2765625,0.4765625,0.5546875,0.3328125,0.5546875,0.5546875,0.5,0.5546875,0.5546875,0.2765625,0.5546875,0.5546875,0.221875,0.240625,0.5,0.221875,0.8328125,0.5546875,0.5546875,0.5546875,0.5546875,0.3328125,0.5,0.2765625,0.5546875,0.5,0.721875,0.5,0.5,0.5,0.3546875,0.259375,0.353125,0.5890625]
  const avg = 0.5279276315789471
  return string
    .split('')
    .map(c => c.charCodeAt(0) < widths.length ? widths[c.charCodeAt(0)] : avg)
    .reduce((cur, acc) => acc + cur) * fontSize
}

function point_on_circle(cx, cy, r, theta) {
  return [cx+Math.cos(theta)*r, cy+Math.sin(theta)*r];
}


// MAIN CLASS
class Alignerator {

  constructor(data) {
    console.log('constructing alignerator');
    console.log(data);
    this.refSequence = data.ref_seq.toUpperCase();
    this.refLen = this.refSequence.length;
    this.refSequenceRC = rc(this.refSequence)
    this.refFeatures = data.ref_features;
    console.log(this.refFeatures);
    this.ref_is_genbank = this.refFeatures.length > 0;
    this.alignments = data.alignments;
    this.base_alignments = [];
    this.w = 800; //window.innerWidth*9/10;
    this.inset_w = 1000;
    this.inset_h = 600;
    this.ref_r = this.w/5;
    this.c = this.w/2;
    this.feature_h = 20;
    this.align_r = 4;
    this.align_r_step = 6;
    this.inset_offset = 5;
    this.inset_top = 20;
    this.inset_al_h = 12;
    this.inset_ref_h = 70;
    this.alignment_page = 0;
    this.svg = d3.select('#pLalignViz').append('svg')
      .attr('id', 'alignment_svg')
      .attr('width', this.w)
      .attr('height', this.w)
      .style('left', 100);

    this.inset_svg = d3.select('#pLalignViz').append('svg')
      .attr('width', this.inset_w)
      .attr('height', this.inset_h)
      .attr('id', 'inset_svg')
      .style('top', this.w);

    this.setup_tooltip();
  }

  setup_tooltip() {
    this.tooltip = d3.select('#pLalignViz').append('div')
      .attr('class', 'WOW_tooltip')
      .html('<h2>yeah</h2><p>uhhuh</p>');
  }
  
  show_tooltip(x, y, text) {
    if (x < 208) {
      this.tooltip.style('left', String(x+8)+'px');
    } else {
      this.tooltip.style('left', String(x-208)+'px');
    }
    this.tooltip.style('display', 'block').html(text);
    this.tooltip.style('top', String(y+8)+'px');
  }
  
  hide_tooltip() {
    this.tooltip.style('display', 'none');
  }

  make_feature_path(d) {
    const start_theta = 2*Math.PI*(d.start/this.refLen);
    const end_theta = 2*Math.PI*(d.end/this.refLen);
    const chevron_theta = Math.min((end_theta-start_theta)/10, Math.PI/40);
    let p = ''
    const arc_outside = 'A ' + String(this.ref_r+this.feature_h/2)+' '+String(this.ref_r+this.feature_h/2)+' 0 ';
    const arc_inside = 'A ' + String(this.ref_r-this.feature_h/2)+' '+String(this.ref_r-this.feature_h/2)+' 0 ';
    if (d.strand == 1) {
      let s1 = point_on_circle(this.c, this.c, this.ref_r-this.feature_h/2, start_theta);
      let s2 = point_on_circle(this.c, this.c, this.ref_r+this.feature_h/2, start_theta);
      let e1 = point_on_circle(this.c, this.c, this.ref_r-this.feature_h/2, end_theta-chevron_theta);
      let e2 = point_on_circle(this.c, this.c, this.ref_r+this.feature_h/2, end_theta-chevron_theta);
      let tip = point_on_circle(this.c, this.c, this.ref_r, end_theta);
      p += 'M'+String(s1[0])+' '+String(s1[1])+' '
      p += 'L'+String(s2[0])+' '+String(s2[1])+' '
      p += arc_outside + '0 1 '+String(e2[0])+' '+String(e2[1])+' '
      p += 'L'+String(tip[0])+' '+String(tip[1])+' '
      p += 'L'+String(e1[0])+' '+String(e1[1])+' '
      p += arc_inside + '0 0 '+String(s1[0])+' '+String(s1[1])+' '
    } else {
      let s1 = point_on_circle(this.c, this.c, this.ref_r-this.feature_h/2, start_theta+chevron_theta);
      let s2 = point_on_circle(this.c, this.c, this.ref_r+this.feature_h/2, start_theta+chevron_theta);
      let e1 = point_on_circle(this.c, this.c, this.ref_r-this.feature_h/2, end_theta);
      let e2 = point_on_circle(this.c, this.c, this.ref_r+this.feature_h/2, end_theta);
      let tip = point_on_circle(this.c, this.c, this.ref_r, start_theta);
      p += 'M'+String(e1[0])+' '+String(e1[1])+' '
      p += 'L'+String(e2[0])+' '+String(e2[1])+' '
      p += arc_outside + '0 0 '+String(s2[0])+' '+String(s2[1])+' '
      p += 'L'+String(tip[0])+' '+String(tip[1])+' '
      p += 'L'+String(s1[0])+' '+String(s1[1])+' '
      p += arc_inside + '0 1 '+String(e1[0])+' '+String(e1[1])+' '
    }
    return p;
  }

  make_feature_text_path(d) {
    const start_theta = 2*Math.PI*(d.start/this.refLen);
    const end_theta = 2*Math.PI*(d.end/this.refLen);
    const chevron_theta = Math.min((end_theta-start_theta)/10, Math.PI/40);
    const arc_inside = 'A ' + String(this.ref_r-this.feature_h/2+3)+' '+String(this.ref_r-this.feature_h/2+3)+' 0 ';
    const s1 = point_on_circle(this.c, this.c, this.ref_r-this.feature_h/2+3, start_theta+chevron_theta);
    const e1 = point_on_circle(this.c, this.c, this.ref_r-this.feature_h/2+3, end_theta-chevron_theta);
    let p = 'M'+String(s1[0])+' '+String(s1[1])+' '
    p += arc_inside + '0 1 '+String(e1[0])+' '+String(e1[1])+' '
    return p;
  }

  make_arc(start_theta, end_theta, r, clockwise_flag, large_arc_flag) {
    const s1 = point_on_circle(this.c, this.c, r, start_theta);
    const e1 = point_on_circle(this.c, this.c, r, end_theta);
    let p = 'M'+String(s1[0])+' '+String(s1[1])+' '
    p += 'A '+ String(r) + ' ' + String(r) + ' 0 ';
    p += String(large_arc_flag) + ' ' + String(clockwise_flag) + ' '; 
    p += String(e1[0])+' '+String(e1[1])+' ';
    return p;
  }

  make_alignment_arc(refpos1, refpos2, r, strand) {
    const alignment_len = Math.abs(refpos2-refpos1);
    const clockwise_flag = strand == 1 ? 1 : 0;
    const large_arc_flag = alignment_len > this.refLen/2 ? 1 : 0;
    const start_theta = 2*Math.PI*(refpos1/this.refLen);
    const end_theta = 2*Math.PI*(refpos2/this.refLen);
    return this.make_arc(start_theta, end_theta, r, clockwise_flag, large_arc_flag);
  }

  setup_drag() {
    const self = this;
    this.drag_r = this.ref_r-50;

    this.drag_el = this.svg.append('g');

    this.drag_circle = this.drag_el.append('circle')
      .attr('cx', self.c)
      .attr('cy', self.c)
      .attr('r', this.drag_r)
      .attr('stroke', '#333')
      .attr('fill', 'none')
      .attr('stroke-width', 30);

    this.drag_arc = this.drag_el.append('path')
      .attr('fill', 'none')
      .attr('stroke', 'red')
      .attr('opacity', 0.8)
      .attr('stroke-width', 30);

    // I realized D3 has an arc generator after I coded it up myself
    // for the alignments...
    const arcGenerator = d3.arc()
      .innerRadius(this.ref_r-50)
      .outerRadius(this.ref_r-50); 

    const dragBehavior = d3.drag()
      .on("start", dragStart)
      .on("drag", drag)
      .on("end", dragEnd);

    this.drag_el.call(dragBehavior);

    let startAngle;
    let clockwise;

    function getAngle(x, y) {
      const at2 = Math.atan2(y, x);
      if (at2 < 0) {
        return (Math.PI*2+at2)
      } else {
        return at2
      }
    }

    function dragStart(e) {
      const [x, y] = [e.x, e.y];
      startAngle = getAngle(x - self.c, y - self.c);
      clockwise = ['not sure yet', true];
    }

    function drag(e) {
      const [x, y] = [e.x, e.y];
      const currentAngle = getAngle(x - self.c, y - self.c);
      if (clockwise[0] === 'not sure yet') {
        // edge case
        let useAngle;
        if (Math.abs(currentAngle-startAngle)>Math.PI) {
          useAngle = (currentAngle > startAngle) ? currentAngle-Math.PI*2 : currentAngle+Math.PI*2;
        } else {
          useAngle = currentAngle;
        }
        if (useAngle > startAngle+Math.PI/10) { // arbitrary condition
          clockwise = ['locked_in', 1];
        } else if (useAngle < startAngle-Math.PI/10) {
          clockwise = ['locked_in', 0];
        } else {
          clockwise = ['not sure yet', useAngle > startAngle ? 1 : 0]
        }
      }
      let angle_difference;
      if ((clockwise[1] == 1) & (currentAngle < startAngle)) {
        angle_difference = Math.abs(Math.PI*2+currentAngle-startAngle);
      } else if ((clockwise[1] == 0) & (currentAngle > startAngle)) {
        angle_difference = Math.abs(Math.PI*-2+currentAngle-startAngle);
      } else {
        angle_difference = Math.abs(currentAngle-startAngle);
      }
      const large_arc = (angle_difference > Math.PI) ? 1 : 0;
      const arcPath = self.make_arc(startAngle, currentAngle, self.drag_r, clockwise[1], large_arc);
      // Update the arc path element (assuming you have one)
      self.drag_arc.attr("d", arcPath);
    }

    function dragEnd(e) {
      const [x, y] = [e.x, e.y];
      const currentAngle = getAngle(x - self.c, y - self.c);
      let start;
      let end;
      if ((clockwise[1] == 1) & (currentAngle < startAngle)) {
        start = startAngle;
        end = currentAngle;
      } else if ((clockwise[1] == 0) & (currentAngle > startAngle)) {
        start = currentAngle;
        end = startAngle;
      } else if (currentAngle > startAngle) {
        start = startAngle;
        end = currentAngle;
      } else {
        start = currentAngle;
        end = startAngle;
      }
      self.set_focal_region(Math.round(self.refLen*start/(2*Math.PI)), 
                            Math.round(self.refLen*end/(2*Math.PI)));
    }
  }

  redraw_inner_arc() {
    const startAngle = 2*Math.PI*this.focal_start/this.refLen;
    const endAngle = 2*Math.PI*this.focal_end/this.refLen;
    const large_arc = (this.domain_len > this.refLen/2) ? 1 : 0;
    const arcPath = this.make_arc(startAngle, endAngle, this.drag_r, 1, large_arc);
    // Update the arc path element (assuming you have one)
    this.drag_arc.attr("d", arcPath);
  }

  set_focal_region(start_base, end_base) {
    this.focal_start = start_base;
    this.focal_end = end_base;
    this.domain_len = (this.focal_start < this.focal_end) ? this.focal_end-this.focal_start : this.refLen-this.focal_start+this.focal_end;
    this.redraw_inner_arc();
    this.draw_focal_region();
  }

  zoom_in() {
    this.focal_start = (this.focal_start + Math.round(this.domain_len/4)) % this.refLen;
    this.focal_end = this.focal_end - Math.round(this.domain_len/4);
    if (this.focal_end < 0) this.focal_end = this.refLen + this.focal_end;
    this.set_focal_region(this.focal_start, this.focal_end);
  }

  zoom_out() {
    const add_len = (this.domain_len < this.refLen/2) ? Math.floor(this.domain_len/2) : Math.floor((this.refLen - this.domain_len) / 2)-1;
    this.focal_start = this.focal_start - add_len;
    if (this.focal_start < 0) this.focal_start = this.refLen + this.focal_start;
    this.focal_end = (this.focal_end + add_len) % this.refLen;
    this.set_focal_region(this.focal_start, this.focal_end);
  }

  setup_focal_region() {
    self = this;
    this.midpoint = this.inset_w/2;

    this.gene_g = this.inset_svg.append('g');

    this.drag_zone = this.inset_svg.append('rect')
      .attr('width', self.inset_w)
      .attr('y', self.inset_ref_h-20)
      .attr('height', 20)
      .attr('fill', '#333');

    this.drag_rect = this.inset_svg.append('rect')
      .attr('fill', 'red')
      .attr('opacity', 0.8)
      .attr('y', self.inset_ref_h-20)
      .attr('height', 20)
      .attr('x', 0)
      .attr('width', 0)

    this.zoom_in_g = this.inset_svg.append('g')
      .attr('class', 'zoom_thing')
      .on('click', (e) => self.zoom_in())
    
    this.zoom_out_g = this.inset_svg.append('g')
      .attr('class', 'zoom_thing')
      .on('click', (e) => self.zoom_out())
    
    this.zoom_in_g.append('circle')
      .attr('cx', self.midpoint+17)
      .attr('cy', 7)
      .attr('r', 6)
      .attr('stroke', base_color)

    this.zoom_in_g.append('line')
      .attr('x1', this.midpoint+17-3)
      .attr('x2', this.midpoint+17+3)
      .attr('y1', 7)
      .attr('y2', 7)
      .attr('stroke', base_color)

    this.zoom_in_g.append('line')
      .attr('x1', this.midpoint+17)
      .attr('x2', this.midpoint+17)
      .attr('y1', 7-3)
      .attr('y2', 7+3)
      .attr('stroke', base_color)

    this.zoom_out_g.append('circle')
      .attr('cx', this.midpoint-17)
      .attr('cy', 7)
      .attr('r', 6)
      .attr('stroke', base_color)

    this.zoom_out_g.append('line')
      .attr('x1', this.midpoint-17-3)
      .attr('x2', this.midpoint-17+3)
      .attr('y1', 7)
      .attr('y2', 7)
      .attr('stroke', base_color)


    const dragAction = d3.drag()
      .on('start', function(e) {
        self.inset_drag_start_mouse = e.x;
      })
      .on('drag', function(e) {
        let start = Math.max(Math.min(e.x, self.inset_drag_start_mouse), 0);
        let end = Math.min(Math.max(e.x, self.inset_drag_start_mouse), self.inset_w);
        self.drag_rect.attr('x', start).attr('width', end-start);
      })
      .on('end', function(e) {
        let start = Math.max(Math.min(e.x, self.inset_drag_start_mouse), 0);
        let end = Math.min(Math.max(e.x, self.inset_drag_start_mouse), self.inset_w);
        self.drag_rect.attr('x', start).attr('width', 0);
        self.set_focal_region(self.focal_start+Math.floor(start*self.domain_len/self.inset_w), self.focal_start+Math.ceil(end*self.domain_len/self.inset_w));
      })
    this.drag_zone.call(dragAction);
  }

  make_gene_display(d) {
    let left;
    let right;
    if (d.indexed_start < this.domain_len) {
      left = Math.round(this.inset_w*d.indexed_start/this.domain_len);
      right = Math.round(this.inset_w*d.indexed_end/this.domain_len);
    } else {
      left = Math.round(this.inset_w*(d.indexed_start-this.refLen)/this.domain_len);
      right = Math.round(this.inset_w*d.indexed_end/this.domain_len);
    }
    const width = right-left;
    const height = 20; //Math.max(Math.min(30, 1000000/this.domain_len), 20);
    const halfHeight = height / 2;
    const chevron_size = (width < 10) ? 0 : Math.min(width/4, 20);
    const top = 20
    const points = d.strand === -1 ? `${left},${top+halfHeight} ${left+chevron_size},${top+height} ${left+width},${top+height} ${left+width},${top} ${left+chevron_size},${top}` : `${right},${top+halfHeight} ${right-chevron_size},${top+height} ${right-width},${top+height} ${right-width},${top} ${right-chevron_size},${top}`;

    const fontsizes = [10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]
    const textBuf = 2.5
    let label = d.label
    let fontsize = fontsizes[0]
    let labelsize = measureText(label, fontsize)
    let labelVisible = (labelsize+2*textBuf+chevron_size < right-left)
    if (labelVisible) {
      for (let f of fontsizes) {
        labelsize = measureText(label, f)
        if (labelsize+2*textBuf+chevron_size < right-left) {
          fontsize = f
        } else {
          break
        }
      }
    }
    const x_pos = d.strand === 1 ? left+textBuf : left+textBuf+chevron_size;
    const y_pos = top+height-textBuf-2;
    const stroke = base_color;
    const strokeWid = 1;
    const chev = `<polygon points="${points}" stroke=${stroke} fill="#333" stroke-width=${strokeWid} />`
    const label_use = labelVisible ? `<text x=${x_pos} y=${y_pos} fill="#FFF">${label}</text>` : '';
    return chev+label_use;
  }

  gene_filter(d) {
    return !((d.indexed_start > this.domain_len) && (d.indexed_end > this.domain_len) && (d.indexed_end > d.indexed_start))
  }

  reindex_genes() {
    for (let g of this.refFeatures) {
      g.indexed_start = (g.start - this.focal_start + this.refLen) % this.refLen;
      g.indexed_end = (g.end - this.focal_start + this.refLen) % this.refLen;
    }
  }

  draw_inset_genes() {
    this.gene_g.remove();
    this.gene_g = this.inset_svg.append('g');
    this.reindex_genes();
    this.gene_g.selectAll('.sbg_gene')
      .data(self.refFeatures.filter((d) => self.gene_filter(d)))
      .enter()
      .append('g')
        .attr('class', 'sgb_gene')
        .attr('opacity', 0.8)
        .html((d) => self.make_gene_display(d))
  }

  draw_focal_region() {
    d3.selectAll('.zoom_thing').style('display', 'block')
    const self = this;
    const base_w = this.inset_w/this.domain_len;
    let bases = [];
    for (let i=0; i<this.domain_len; i++) {
      bases.push((this.focal_start+i) % this.refLen);
    }

    this.inset_svg.selectAll('.ref_base').remove();
    this.inset_svg.selectAll('.ref_base')
      .data(bases)
      .enter()
      .append('g')
        .attr('class', 'ref_base')

    this.inset_svg.selectAll('.ref_base')
      .append('rect')
        .attr('x', (d, i) => i*base_w)
        .attr('y', this.inset_top+this.inset_ref_h)
        .attr('width', base_w)
        .attr('height', this.inset_al_h)
        .attr('opacity', 0.7)
        .attr('fill', (d) => colors[self.refSequence[d]]);
    
    if (base_w > 10) {
      this.inset_svg.selectAll('.ref_base')
        .append('text')
          .attr('x', (d, i) => (i+0.5)*base_w)
          .attr('y', this.inset_top+this.inset_ref_h+this.inset_al_h-2)
          .attr('text-anchor', 'middle')
          .attr('fill', 'white')
          .attr('font-size', this.inset_al_h)
          .html((d) => self.refSequence[d]);
    }

    this.draw_inset_genes();
    this.inset_svg.selectAll('.inset_alignment').remove();
    for (let i=0; i<this.base_alignments.length; i++) {
      this.draw_inset_alignment(bases, base_w, i)
    }
  }

  parse_tooltip_text(t) {
    let s = ''
    for (let ab of t) {
      if (ab[0] == 'insertion') {
        s += '<h2>' + String(ab[2].length)+'bp insertion</h2>' + ab[2].slice(0, 10) + ((ab[2].length > 10) ? '...' : '')
      } else {
        s += '<p>'+String(ab[0])+' : '+String(ab[1])+'</p>';
      }
    }
    return s;
  }

  draw_inset_alignment(bases, base_w, al_index) {
    const self = this;
    let al_el = this.inset_svg.append('g').attr('class', 'inset_alignment').style('white-space', 'pre-wrap');
    for (let i=0; i<bases.length; i++) {
  
      let al = this.base_alignments[al_index][bases[i]];
      let g = al_el.append('g')
        .on('mouseover', (e) => self.show_tooltip(e.x, e.y, self.parse_tooltip_text(al)))
        .on('mouseout', (e) => self.hide_tooltip())
        .on('click', function(e) {
          console.log(self.alignments[al_index]);
        });
      /*
      if (al.filter((ab) => ab[0]!='insertion').length > 1) { // overlap
        let consensus = Array(...new Set(al.map((ab) => ab[1])));
        if (consensus.length == 1) {
          al = [['overlap', consensus[0]]];
        } else {
          al = al.concat(al.map((ab) => ['overlap', ab[1]]));
        }
      }
      */
      let al_height_use = (this.inset_al_h / al.length) ;
      let ab_index = 0
      for (let ab of al) {
        let top = this.inset_top+this.inset_ref_h+(al_index+2)*(this.inset_al_h+this.inset_offset)+ab_index*(al_height_use)
        ab_index += 1
        g.append('rect')
          .attr('x', i*base_w)
          .attr('y', top)
          .attr('width', base_w)
          .attr('height', al_height_use) //(ab[0]=='insertion') ? this.inset_al_h+5 : this.inset_al_h)
          .attr('fill', (ab[0]=='mismatch') ? colors[ab[1]] : colors[ab[0]])
          .attr('title', ab.join(' '));

        if (base_w > 10) {
          g.append('text')
              .attr('x', (i+0.5)*base_w)
              .attr('y', top+al_height_use-2)
              .attr('text-anchor', 'middle')
              .attr('fill', 'white')
              .attr('font-size', al_height_use)
              .html((d) => ab[1]);
        }
      }
    }
  }

  draw_the_ref() {
    const self = this;
    this.setup_drag();

    this.svg.append('circle')
      .attr('cx', self.c)
      .attr('cy', self.c)
      .attr('r', this.ref_r)
      .attr('stroke', '#CCC')
      .attr('fill', 'none')
      .attr('stroke-width', 4);

    if (this.ref_is_genbank) {
      this.svg.selectAll('.genbank_feature')
        .data(self.refFeatures.filter((d) => (d.end-d.start)>50))
        .enter()
        .append('g')
          .attr('class', 'genbank_feature');

      this.svg.selectAll('.genbank_feature').append('path')
        .attr('class', 'genbank_feature_path')
        .attr('stroke', '#FFF')
        .attr('fill', '#999')
        .attr('d', (d) => self.make_feature_path(d));

      this.svg.selectAll('.genbank_feature').append('path')
        .attr('class', 'genbank_feature_text_path')
        .attr('fill', 'none')
        .attr('id', (d, i) => 'text_path'+String(i))
        .attr('d', (d) => self.make_feature_text_path(d))

      this.svg.selectAll('.genbank_feature').append('text').append('textPath')
        .attr('href', (d, i) => '#text_path'+String(i))
        .html((d) => d.label)
    }
    this.setup_focal_region();
  }

  increment_alignment_page(step) {
    this.alignment_page += step;
    console.log(this.alignment_page, this.alignments.length, Math.floor(this.alignments.length/alignment_num));
    const on_last_page = (this.alignment_page >= Math.floor((this.alignments.length-1)/alignment_num));
    this.next_button.style('display', on_last_page ? 'none': 'block');
    this.prev_button.style('display', this.alignment_page == 0 ? 'none' : 'block');
    this.draw_alignments();
    this.draw_focal_region();
  }

  draw_al_page_buttons() {
    const self = this;
    this.next_button = this.svg.append('text')
      .attr('class', 'next_prev_button')
      .attr('x', self.w-20)
      .attr('y', self.w-20)
      .attr('text-anchor', 'end')
      .style('display', self.alignments.length <= alignment_num ? 'none' : 'block')
      .attr('fill', colors['button_text'])
      .html('next '+String(alignment_num) + ' ->')
      .on('click', () => self.increment_alignment_page(1));
    this.prev_button = this.svg.append('text')
      .attr('class', 'next_prev_button')
      .attr('x', 20)
      .attr('y', self.w-20)
      .attr('text-anchor', 'start')
      .style('display', 'none')
      .attr('fill', colors['button_text'])
      .html('<- previous')
      .on('click', () => self.increment_alignment_page(-1));
  }

  refpos(strand, i) {
    while (i > (this.refLen-1)) {
      i = i - this.refLen;
    }
    const pos = (i < 0) ? this.refLen+i : i;
    const res = (strand == 1) ? pos : this.refLen-pos-1;
    return res;
  }

  qrc(strand, s) {
    return (strand == 1) ? s : rc(s);
  }

  draw_alignment_old(alignment, alignment_count) {
    // deletions are black arcs
    // insertions are blue lines with arcs representing insertion size
    // mismatches are red arcs
    const self = this;
    const strand = alignment.strand;
    const query = alignment.query_seq;
    
    // flipping things around if it is on the other strand
    // (not the simplest way to do this, but it works)
    const use_seq = strand == 1 ? this.refSequence+this.refSequence : this.refSequenceRC+this.refSequenceRC;
    const ref_start = strand == 1 ? alignment.ref_start : this.refLen*2-alignment.ref_end;
    const use_cigar = strand == 1 ? alignment.cigar : alignment.cigar.slice().reverse();

    let query_position = alignment.query_start;
    let ref_position = ref_start;
    let mismatches = [];
    let deletions = [];
    let insertions = [];
    let alignment_by_ref_pos = [];
    for (let i=0; i<this.refLen; i++) {
      alignment_by_ref_pos.push([]);
    }
    if (query_position > 0) {
      alignment_by_ref_pos[this.refpos(strand, ref_position-1)].push(['pre-alignment', strand == 1 ? '<' : '>', this.qrc(strand, query.slice(0, query_position))]);
    }
    for (let cig of use_cigar) {
      let bp = cig[0];
      let letter = "MIDNSHP=XB"[cig[1]];
      if (letter == 'S') {
        alignment_by_ref_pos[this.refpos(strand, ref_position-1)].push(['pre-alignment', strand == 1 ? '>' : '<', this.qrc(strand, query.slice(query_position, query_position+bp))]);
        query_position += bp;
      } else if (letter == 'M') {
        for (let i=0; i<bp; i++) {
          if (query[query_position+i]!=use_seq[ref_position+i]) {
            mismatches.push(ref_position+i);
            alignment_by_ref_pos[this.refpos(strand, ref_position+i)].push(['mismatch', this.qrc(strand, query[query_position+i])]);
          } else {
            alignment_by_ref_pos[this.refpos(strand, ref_position+i)].push(['match', this.qrc(strand, query[query_position+i])]);
          }
        }
        query_position += bp;
        ref_position += bp;
      } else if (letter == 'D') {
        for (let i=0; i<bp; i++) {
          alignment_by_ref_pos[this.refpos(strand, ref_position+i)].push(['deletion', '-']);
        }
        ref_position += bp;
        deletions.push([ref_position-bp, ref_position]);
      } else if (letter == 'I') {
        alignment_by_ref_pos[this.refpos(strand, ref_position)].push(['insertion', '', this.qrc(strand, query.slice(query_position, query_position+bp))]);
        query_position += bp;
        insertions.push([ref_position, bp]);
      }
    }
    if (query_position < query.length-1) {
      console.log('Q', query.length, query_position);
      alignment_by_ref_pos[this.refpos(strand, ref_position)].push(['post-alignment', strand == 1 ? '>' : '<', this.qrc(strand, query.slice(query_position, query.length))]);
    }
    console.log(ref_start, ref_position, ref_start - ref_position, strand, mismatches.length);
    //console.log(mismatches, deletions, insertions);
    let al_g = this.svg.append('g')
      .attr('class', 'alignment_group');

    let al_color = colors['al'];
    if (Math.abs(ref_start-ref_position) > this.refLen) {
      console.log('alignment longer than reference...');
      // dealing with cases where the alignment is longer than the ref...
      al_g.append('path')
        .attr('class', 'alignment_path')
        .attr('fill', 'none')
        .attr('stroke', colors['al'])
        .attr('stroke-width', this.align_r/2)
        .attr('d', self.make_alignment_arc(ref_start, ref_start+(this.refLen-1)*strand, this.ref_r+this.feature_h+alignment_count*this.align_r_step, strand));

      ref_position = ref_position-this.refLen;
      al_color = colors['overlap'];
    }

    al_g.append('path')
      .attr('class', 'alignment_path')
      .attr('fill', 'none')
      .attr('stroke', al_color)
      .attr('stroke-width', this.align_r/2)
      .attr('d', self.make_alignment_arc(ref_start*strand, ref_position*strand, this.ref_r+this.feature_h+alignment_count*this.align_r_step, strand));

    al_g.selectAll('.mismatch_marks')
      .data(mismatches)
      .enter()
      .append('path')
        .attr('class', 'mismatch_marks')
        .attr('fill', 'none')
        .attr('stroke', colors['X'])
        .attr('stroke-width', this.align_r/2)
        .attr('d', (d) => self.make_alignment_arc(d*strand-0.5, d*strand+0.5, this.ref_r+this.feature_h+alignment_count*this.align_r_step, strand));

    al_g.selectAll('.deletion_arcs')
      .data(deletions)
      .enter()
      .append('path')
        .attr('class', 'deletion_arcs')
        .attr('fill', 'none')
        .attr('stroke', colors['D'])
        .attr('stroke-width', this.align_r/2)
        .attr('d', (d) => self.make_alignment_arc(d[0]*strand-0.5, d[1]*strand+0.5, this.ref_r+this.feature_h+alignment_count*this.align_r_step, strand));
    
    al_g.selectAll('.insertion_arcs')
      .data(insertions)
      .enter()
      .append('path')
        .attr('class', 'insertion_arcs')
        .attr('fill', 'none')
        .attr('stroke', colors['I'])
        .attr('stroke-width', 3*this.align_r/4) //(d) => (d[1]-d[0]) > 5 ? 3*this.align_r/4 : this.align_r/3)
        .attr('d', (d) => self.make_alignment_arc((d[0]-d[1]/2)*strand-0.5, (d[0]+d[1]/2)*strand+0.5, this.ref_r+this.feature_h+alignment_count*this.align_r_step+this.align_r/2, strand));

    return alignment_by_ref_pos
  }

  draw_alignments() {
    d3.selectAll('.alignment_group').remove();
    this.base_alignments = [];
    for (let i=alignment_num*this.alignment_page; 
         i<Math.min(alignment_num*(1+this.alignment_page), this.alignments.length); i++) {
      this.base_alignments.push(this.draw_alignment_old(this.alignments[i], i-alignment_num*this.alignment_page));
    }
    console.log(this.base_alignments);
  }

}

function goforit() {
  const al = new Alignerator(pLalignData);
  al.draw_the_ref();
  al.draw_alignments();
  al.draw_al_page_buttons();
}

goforit();
