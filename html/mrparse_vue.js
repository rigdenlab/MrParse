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
  	<div v-if="ss_pred || ss_pred" id='classification1'>
	    <h2>Classification</h2>
	    <!_- <div id='classification'></div> -->
	    <pfam-region :id="'ss_pred'" :region="ss_pred"/>
	    <pfam-region :id="'classification'" :region="classification"/>
	</div>
	<div v-else id='classification1'>
	    <h4>### Classification step omitted ###</h4>
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
    <th>name</th>
    <th>Resolution</th>
    <th>Space Group</th>
    <th>Has NCS?</th>
    <th>Has Twinning?</th>
    <th>Has Anisotropy?</th>
    <th>File Path</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td>{{ hklinfo.name }}</td>
    <td>{{ hklinfo.resolution | decimalPlaces }}</td>
    <td>{{ hklinfo.space_group }}</td>
    <td>{{ hklinfo.has_ncs }}</td>
    <td>{{ hklinfo.has_twinning }}</td>
    <td>{{ hklinfo.has_anisotropy }}</td>
    <td>{{ hklinfo.hklin }}</td>
  </tr>
</tbody>
</table>
</div>`
});


Vue.component('homolog-table', {
  data: function() {
    return {
      homologs: this.$root.homologs,
      /* Need to write data structure to link columns, titles and tooltips */
      columns: ['name', 'pdb_id', 'chain_id', 'region_id', 'range', 'length', 'ellg', 'molecular_weight', 'rmsd', 'seq_ident'],
      columnTitles: {
        'name': 'Name of the homolog (<PDB>_<CHAIN_ID>_<NUMBER>',
        'pdb_id': 'PDB code of homolog',
        'chain_id': 'Chain ID in PDB',
        'region_id': 'Number of the region',
        'range': 'Start - stop coordinates of the homolog',
        'length': 'Length of the homolog in residues',
        'ellg': 'Computed Log Likelihood Gain',
        'molecular_weight': 'Molecular Weight in Daltons',
        'rmsd': 'RMSD from template',
        'seq_ident': 'Sequence Identity to template'
      },
      sortKey: 'domain',
      order: 'asc',
  	  attrs_to_column_title: {'name' : 'Name',
  	                          'pdb_id' : 'PDB',
							  'chain_id' : 'Chain',
                              'region_id' : 'Region',
                              'range' : 'Range',
                              'length' : 'Length',
                              'ellg' : 'eLLG',
                              'molecular_weight' : 'Mol. Wt.',
                              'rmsd' : 'RMSD',
                              'seq_ident' : 'Seq. Ident.'},
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
    getTitle: function(column) {
      return this.columnTitles[column];
    }
  },
  mounted: function() {
    this.sortBy('region')
  },
  template: `
  <div class="homolog-table">
  <table id="homologs">
      <thead>
        <tr>
          <th v-for="column in columns">
            <a href="#" @click="sortBy(column)" v-bind:title="getTitle(column)">
              {{ attrs_to_column_title[column] }}
            </a>
          </th>
        </tr>
      </thead>
      <tbody>
        <tr v-for="homolog in homologs">
          <td><a v-bind:href="homolog.pdb_file">{{ homolog.name }}</a></td>
          <td><a v-bind:href="homolog.pdb_url" target="_blank">{{ homolog.pdb_id }}</a></td>
          <td>{{ homolog.chain_id }}</td>
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
