/*
data attribute creates object that is bound to this

computed - computed property but with caching
method - computed property - invoked each time

v-model - links inputs to vue js data


https://codepen.io/pespantelis/pen/ojwgPB
https://www.raymondcamden.com/2018/02/08/building-table-sorting-and-pagination-in-vuejs

*/

function add_graphic(region, parent, residueWidth = 1.0) {
  let pg = new PfamGraphic();
  pg.setImageParams({
    residueWidth: residueWidth,
    xscale: 1.0,
    yscale: 1.0
  });
  pg.setParent(parent);
  pg.setSequence(region);
  pg.render();
}


const EventBus = new Vue();


Vue.filter("decimalPlaces", (value, num = 2) => {
    if (value == null) {
        return "N/A";
    } else {
        return value.toFixed(num);
    }
  });


Vue.component('pfam-graphics', {
  data: function() {
    return {
      homologs: this.$root.homologs,
      ss_pred: this.$root.ss_pred,
      classification: this.$root.classification
    }
  },
  created: function() {
    EventBus.$on("sortedData", sortedData => {
      this.homologs = sortedData;
    })
  },
  template: `
  <div class="pfam-graphics" ref=pfamgraphics>
  	<div v-if="ss_pred || ss_pred" id='classification'>
	    <h2>Sequence Based Predictions</h2>
	    <pfam-region :id="'ss_pred'" :region="ss_pred"/>
	    <pfam-region :id="'classification'" :region="classification"/>
	</div>
	<div v-else id='classification'>
	    <h4>### Sequence Based Prediction step omitted ###</h4>
	</div>
	<h2>Regions</h2>
    <pfam-region v-for="homolog in homologs" :key="homolog.name" :id="homolog.name" :region="homolog._pfam_json"/>
  </div>
  `
});


Vue.component('pfam-region', {
  props: {
    region: Object
  },
  mounted: function() {
    let residueWidth = Math.max(1.0, this.$parent.$refs.pfamgraphics.clientWidth / this.region.length);
    add_graphic(this.region, this.$refs.cdiv, residueWidth);
  },
  template: '<div id=id ref=cdiv></div>'
});

Vue.component('hkl-info-table', {
  data: function() {
    return {
      hklinfo: this.$root.hklinfo
    }
  },
  template: `<div v-if="hklinfo" id="hkl_info">
<h2>HKL Info</h2>
<table>
<thead>
  <tr style="text-align: right;">
    <th title='Name of, and link to, the file crystallographic data file.'>name</th>
    <th title='Highest resolution of the crystallographic data'>Resolution</th>
    <th title='The space group of the crystallographic data'>Space Group</th>
    <th title='Indicates the presences of Non-Crystallographic Symmetry (as calcualted by CTRUNCATE)'>Has NCS?</th>
    <th title='Indicates the presences of Twinning (as calcualted by CTRUNCATE)'>Has Twinning?</th>
    <th title='Indicates the presences of Anisotropy (as calcualted by CTRUNCATE)'>Has Anisotropy?</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td><a v-bind:href="hklinfo.hklin">{{ hklinfo.name }}</a></td>
    <td>{{ hklinfo.resolution | decimalPlaces }}</td>
    <td>{{ hklinfo.space_group }}</td>
    <td>{{ hklinfo.has_ncs }}</td>
    <td>{{ hklinfo.has_twinning }}</td>
    <td>{{ hklinfo.has_anisotropy }}</td>
  </tr>
</tbody>
</table>
</div>`
});


Vue.component('homolog-table', {
  data: function() {
    return {
      homologs: this.$root.homologs,
      sortKey: 'domain',
      order: 'asc',
  	  columns: [ { 'attr':  'name',
  	               'title': 'Name',
  	               'popup': 'Name of the homolog (<PDB>_<CHAIN_ID>_<NUMBER>'},
  	             { 'attr': 'pdb_id',
  	              'title': 'PDB',
  	              'popup': 'PDB code of homolog'},   
  	             { 'attr': 'region_id',
  	               'title': 'Region',
  	               'popup': 'Number of the region'},   
  	             { 'attr': 'range',
  	               'title': 'Range',
  	               'popup': 'Start - stop coordinates of the homolog'},   
  	             { 'attr': 'length',
  	               'title': 'Length',
  	                'popup': 'Length of the homolog in residues'},   
  	             { 'attr': 'ellg',
  	               'title': 'eLLG',
  	               'popup': 'Computed Log Likelihood Gain'},   
  	             { 'attr': 'molecular_weight',
  	               'title': 'Mol. Wt.',
  	               'popup': 'Molecular Weight in Daltons'},   
  	             { 'attr': 'rmsd',
  	               'title': 'eRMSD',
  	               'popup': 'Estimated RMSD from template'},   
  	             { 'attr': 'seq_ident',
  	               'title': 'Seq. Ident.',
  	               'popup': 'Sequence Identity to template'}],
    }
  },
  methods: {
    sortBy: function(sortKey) {
      if (this.sortKey == sortKey) {
        if (this.order == 'asc') {
          this.order = 'desc';
        } else {
          this.order = 'asc';
        }
      }
      this.sortKey = sortKey;
      this.homologs = _.orderBy(this.homologs, this.sortKey, this.order);
      EventBus.$emit("sortedData", this.homologs);
      return
    },
  },
  mounted: function() {
    this.sortBy('region_id')
  },
  template: `
  <div class="homolog-table">
  <table id="homologs">
      <thead>
        <tr>
          <th v-for="column in columns">
            <a href="#" @click="sortBy(column.attr)" v-bind:title="column.popup">
              {{ column.title }}
            </a>
          </th>
        </tr>
      </thead>
      <tbody>
        <tr v-for="homolog in homologs">
          <td><a v-bind:href="homolog.pdb_file">{{ homolog.name }}</a></td>
          <td><a v-bind:href="homolog.pdb_url" target="_blank">{{ homolog.pdb_id }}</a></td>
          <td>{{ homolog.region_id }}</td>
          <td>{{ homolog.range }}</td>
          <td>{{ homolog.length }}</td>
          <td>{{ homolog.ellg }}</td>
          <td>{{ homolog.molecular_weight | decimalPlaces }}</td>
          <td>{{ homolog.rmsd }}</td>
          <td>{{ homolog.seq_ident }}</td>
        </tr>
      </tbody>
    </table>
    </div>
    `
})


new Vue({
  el: '#app',
  /* Check hkl_info undefined */
  data: {
    homologs: mrparse_data.pfam.homologs,
    ss_pred: mrparse_data.pfam.ss_pred,
    classification: mrparse_data.pfam.classification,
    hklinfo: mrparse_data.hkl_info
  },
})
