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

/* EventBus is used to pass changes between components */
const EventBus1 = new Vue();
const EventBus2 = new Vue();


Vue.filter("decimalPlaces", (value, num = 2) => {
    if (value == null) {
        return "N/A";
    } else {
        return value.toFixed(num);
    }
});


Vue.component('pfam-graphics', {
    data: function () {
        return {
            homologs: this.$root.homologs,
            ss_pred: this.$root.ss_pred,
            classification: this.$root.classification
        }
    },
    /* When the homolog table has been sorted, the homologs are put on the EventBus so we set this components
       homologs to be the sorted set from the EventBus */
    created: function () {
        EventBus1.$on("sortedData1", sortedData => {
            this.homologs = sortedData;
        })
    },
    template: `
  <div class="pfam-graphics" ref=pfamgraphics>
	<h2 style="font-size:15px;color:#3b3b3dff;font-weight:normal;">Visualisation of Regions</h2>
    <pfam-region v-for="homolog in homologs" :key="homolog.name" :id="homolog.name" :region="homolog._pfam_json"/>
    <div v-if="ss_pred || classification" id    ='classification'>
	    <h2 style="font-size:15px;color:#3b3b3dff;font-weight:normal;">Sequence Based Predictions</h2>
	    <pfam-region :id="'ss_pred'" :region="ss_pred"/>
	    <pfam-region :id="'classification'" :region="classification"/>
	</div>
	<div v-else id='classification'>
	    <h3>** Sequence Based Prediction step was skipped: append <tt>--do_classify</tt> argument to run **</h3>
	</div>
  </div>
  `
});

Vue.component('af-model-pfam-graphics', {
    data: function () {
        return {
            af_models: this.$root.af_models,
        }
    },
    /* When the model table has been sorted, the models are put on the EventBus so we set this components
       models to be the sorted set from the EventBus */
    created: function () {
        EventBus2.$on("sortedData2", sortedData => {
            this.af_models = sortedData;
        })
    },
    template: `
  <div class="af-model-pfam-graphics" ref=afmodelpfamgraphics>
	<h2 style="font-size:15px;color:#3b3b3dff;font-weight:normal;">Visualisation of Regions</h2>
    <af-model-pfam-region v-for="af_model in af_models" :key="af_model.name" :id="af_model.name" :region="af_model._pfam_json"/>
  </div>
  `
});

Vue.component('esm-model-pfam-graphics', {
    data: function () {
        return {
            esm_models: this.$root.esm_models,
        }
    },
    /* When the model table has been sorted, the models are put on the EventBus so we set this components
       models to be the sorted set from the EventBus */
    created: function () {
        EventBus3.$on("sortedData3", sortedData => {
            this.esm_models = sortedData;
        })
    },
    template: `
  <div class="esm-model-pfam-graphics" ref=esmmodelpfamgraphics>
	<h2 style="font-size:15px;color:#3b3b3dff;font-weight:normal;">Visualisation of Regions</h2>
    <esm-model-pfam-region v-for="esm_model in esm_models" :key="esm_model.name" :id="esm_model.name" :region="esm_model._pfam_json"/>
  </div>
  `
});

Vue.component('pfam-region', {
    props: {
        region: Object
    },
    mounted: function () {
        let residueWidth = Math.max(1.0, this.$parent.$refs.pfamgraphics.clientWidth / this.region.length);
        add_graphic(this.region, this.$refs.cdiv1, residueWidth);
    },
    template: '<div id=id ref=cdiv1></div>'
});

Vue.component('af-model-pfam-region', {
    props: {
        region: Object
    },
    mounted: function () {
        let residueWidth = Math.max(1.0, this.$parent.$refs.afmodelpfamgraphics.clientWidth / this.region.length);
        add_graphic(this.region, this.$refs.cdiv2, residueWidth);
    },
    template: '<div id=id ref=cdiv2></div>'
});

Vue.component('esm-model-pfam-region', {
    props: {
        region: Object
    },
    mounted: function () {
        let residueWidth = Math.max(1.0, this.$parent.$refs.esmmodelpfamgraphics.clientWidth / this.region.length);
        add_graphic(this.region, this.$refs.cdiv3, residueWidth);
    },
    template: '<div id=id ref=cdiv3></div>'
});

Vue.component('hkl-info-table', {
    data: function () {
        return {
            hklinfo: this.$root.hklinfo
        }
    },
    template: `<div v-if="hklinfo" id="hkl_info">
<h2>HKL Info</h2>
<div class="hkl-table">
<table>
<thead>
  <tr style="text-align: right;">
    <th title='Name of, and link to, the file crystallographic data file.'>Name</th>
    <th title='Highest resolution of the crystallographic data'>Resolution</th>
    <th title='The space group of the crystallographic data'>Space Group</th>
    <th title='Indicates the presences of Non-Crystallographic Symmetry (as calculated by CTRUNCATE)'>Has NCS?</th>
    <th title='Indicates the presences of Twinning (as calculated by CTRUNCATE)'>Has Twinning?</th>
    <th title='Indicates the presences of Anisotropy (as calculated by CTRUNCATE)'>Has Anisotropy?</th>
  </tr>
</thead>
<tbody>
  <tr>
    <td>
    <a v-if="hklinfo.hklin" :href="hklinfo.hklin">{{ hklinfo.name }}</a>
    <a v-else>{{ hklinfo.name }}</a>
    </td>
    <td>{{ hklinfo.resolution | decimalPlaces }}</td>
    <td>{{ hklinfo.space_group }}</td>
    <td>{{ hklinfo.has_ncs }}</td>
    <td>{{ hklinfo.has_twinning }}</td>
    <td>{{ hklinfo.has_anisotropy }}</td>
  </tr>
</tbody>
</table>
</div>
</div>`
});


Vue.component('homolog-table', {
    data: function () {
        return {
            homologs: this.$root.homologs,
            sortKey: 'domain',
            order: 'desc',
            columns: [{
                'attr': 'name',
                'title': 'Name',
                'popup': 'Name of the homolog (<PDB>_<CHAIN_ID>_<NUMBER>'
            },
                {
                    'attr': 'pdb_id',
                    'title': 'PDB',
                    'popup': 'PDB code of homolog'
                },
                {
                    'attr': 'resolution',
                    'title': 'Resolution',
                    'popup': 'Highest resolution of Xtal data for PDB'
                },
                {
                    'attr': 'region_id',
                    'title': 'Region',
                    'popup': 'Number of the region'
                },
                {
                    'attr': 'range',
                    'title': 'Range',
                    'popup': 'Start - stop coordinates of the homolog'
                },
                {
                    'attr': 'length',
                    'title': 'Length',
                    'popup': 'Length of the homolog in residues'
                },
                {
                    'attr': 'ellg',
                    'title': 'eLLG',
                    'popup': 'Computed Log Likelihood Gain'
                },
                {
                    'attr': 'molecular_weight',
                    'title': 'Mol. Wt.',
                    'popup': 'Molecular Weight in Daltons'
                },
                {
                    'attr': 'rmsd',
                    'title': 'eRMSD',
                    'popup': 'Estimated RMSD from template'
                },
                {
                    'attr': 'seq_ident',
                    'title': 'Seq. Ident.',
                    'popup': 'Sequence Identity to template'
                }],
        }
    },
    methods: {
        sortBy: function (sortKey) {
            if (this.sortKey == sortKey) {
                if (this.order == 'asc') {
                    this.order = 'desc';
                } else {
                    this.order = 'asc';
                }
            }
            this.sortKey = sortKey;
            /* sorting is done using _.orderBy from the loadsh libary */
            this.homologs = _.orderBy(this.homologs, this.sortKey, this.order);
            /* We've sorted the homlogs, so put them on the EventBus so that the graphics will be updated */
            EventBus1.$emit("sortedData1", this.homologs);
            return
        },
    },
    mounted: function () {
        this.sortBy('seq_ident')
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
          <td>
          <a v-if="homolog.pdb_file" :href="homolog.pdb_file">{{ homolog.name }}</a>
          <a v-else>{{ homolog.name }}</a>
          </td>
          <td><a v-bind:href="homolog.pdb_url" target="_blank">{{ homolog.pdb_id }}</a></td>
          <td>{{ homolog.resolution  | decimalPlaces }}</td>
          <td>{{ homolog.region_id }}</td>
          <td>{{ homolog.range }}</td>
          <td>{{ homolog.length }}</td>
          <td>{{ homolog.ellg }}</td>
          <td>{{ homolog.molecular_weight }}</td>
          <td>{{ homolog.rmsd }}</td>
          <td>{{ homolog.seq_ident }}</td>
        </tr>
      </tbody>
    </table>
    </div>
    `
})

Vue.component('af-model-table', {
    data: function () {
        return {
            af_models: this.$root.af_models,
            sortKey: 'domain',
            order: 'desc',
            columns: [{
                'attr': 'name',
                'title': 'Name',
                'popup': 'Name of the model'
            },
                {
                    'attr': 'model_id',
                    'title': 'model',
                    'popup': 'Code of model'
                },
                {
                    'attr': 'date_made',
                    'title': 'Date Made',
                    'popup': 'The date the model was made (or released)'
                },
                {
                    'attr': 'region_id',
                    'title': 'Region',
                    'popup': 'Number of the region'
                },
                {
                    'attr': 'range',
                    'title': 'Range',
                    'popup': 'Start - stop coordinates of the model'
                },
                {
                    'attr': 'length',
                    'title': 'Length',
                    'popup': 'Length of the model in residues'
                },
                {
                    'attr': 'avg_plddt',
                    'title': 'Avg. pLDDT',
                    'popup': 'The average pLDDT score in the model'
                },
                {
                    'attr': 'h_score',
                    'title': 'H-score',
                    'popup': 'H-index style quality score measuring percentage of model vs pLDDT'
                },
                {
                    'attr': 'seq_ident',
                    'title': 'Seq. Ident.',
                    'popup': 'Sequence Identity to template'
                }],
        }
    },

    methods: {
        sortBy: function (sortKey) {
            if (this.sortKey == sortKey) {
                if (this.order == 'asc') {
                    this.order = 'desc';
                } else {
                    this.order = 'asc';
                }
            }
            this.sortKey = sortKey;
            /* sorting is done using _.orderBy from the loadsh libary */
            this.af_models = _.orderBy(this.af_models, this.sortKey, this.order);
            /* We've sorted the models, so put them on the EventBus so that the graphics will be updated */
            EventBus2.$emit("sortedData2", this.af_models);
            return
        },
    },
    mounted: function () {
        this.sortBy('seq_ident')
    },
    template: `
  <div class="af-model-table">
  <table id="af-models">
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
        <tr v-for="afmodel in af_models">
          <td>
          <a v-if="afmodel.pdb_file" :href="afmodel.pdb_file">{{ afmodel.name }}</a>
          <a v-else>{{ afmodel.name }}</a>
          </td>
          <td><a v-bind:href="afmodel.model_url" target="_blank">{{ afmodel.model_id }}</a></td>
          <td>{{ afmodel.date_made }}</td>
          <td>{{ afmodel.region_id }}</td>
          <td>{{ afmodel.range }}</td>
          <td>{{ afmodel.length }}</td>
          <td>{{ afmodel.avg_plddt | decimalPlaces }}</td>
          <td>{{ afmodel.h_score }}</td>
          <td>{{ afmodel.seq_ident }}</td>
        </tr>
      </tbody>
    </table>
    </div>
    `
})

Vue.component('esm-model-table', {
    data: function () {
        return {
            esm_models: this.$root.esm_models,
            sortKey: 'domain',
            order: 'desc',
            columns: [{
                'attr': 'name',
                'title': 'Name',
                'popup': 'Name of the model'
            },
                {
                    'attr': 'model_id',
                    'title': 'model',
                    'popup': 'Code of model'
                },
                {
                    'attr': 'date_made',
                    'title': 'Date Made',
                    'popup': 'The date the model was made (or released)'
                },
                {
                    'attr': 'region_id',
                    'title': 'Region',
                    'popup': 'Number of the region'
                },
                {
                    'attr': 'range',
                    'title': 'Range',
                    'popup': 'Start - stop coordinates of the model'
                },
                {
                    'attr': 'length',
                    'title': 'Length',
                    'popup': 'Length of the model in residues'
                },
                {
                    'attr': 'avg_plddt',
                    'title': 'Avg. pLDDT',
                    'popup': 'The average pLDDT score in the model'
                },
                {
                    'attr': 'h_score',
                    'title': 'H-score',
                    'popup': 'H-index style quality score measuring percentage of model vs pLDDT'
                },
                {
                    'attr': 'seq_ident',
                    'title': 'Seq. Ident.',
                    'popup': 'Sequence Identity to template'
                }],
        }
    },

    methods: {
        sortBy: function (sortKey) {
            if (this.sortKey == sortKey) {
                if (this.order == 'asc') {
                    this.order = 'desc';
                } else {
                    this.order = 'asc';
                }
            }
            this.sortKey = sortKey;
            /* sorting is done using _.orderBy from the loadsh libary */
            this.esm_models = _.orderBy(this.esm_models, this.sortKey, this.order);
            /* We've sorted the models, so put them on the EventBus so that the graphics will be updated */
            EventBus3.$emit("sortedData3", this.esm_models);
            return
        },
    },
    mounted: function () {
        this.sortBy('seq_ident')
    },
    template: `
  <div class="esm-model-table">
  <table id="esm-models">
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
        <tr v-for="esmmodel in esm_models">
          <td>
          <a v-if="esmmodel.pdb_file" :href="esmmodel.pdb_file">{{ esmmodel.name }}</a>
          <a v-else>{{ esmmodel.name }}</a>
          </td>
          <td><a v-bind:href="esmmodel.model_url" target="_blank">{{ esmmodel.model_id }}</a></td>
          <td>{{ esmmodel.date_made }}</td>
          <td>{{ esmmodel.region_id }}</td>
          <td>{{ esmmodel.range }}</td>
          <td>{{ esmmodel.length }}</td>
          <td>{{ esmmodel.avg_plddt | decimalPlaces }}</td>
          <td>{{ esmmodel.h_score }}</td>
          <td>{{ esmmodel.seq_ident }}</td>
        </tr>
      </tbody>
    </table>
    </div>
    `
})

new Vue({
    el: '#app',
    data: {
        homologs: mrparse_data.pfam.homologs,
        ss_pred: mrparse_data.pfam.ss_pred,
        classification: mrparse_data.pfam.classification,
        hklinfo: mrparse_data.hkl_info,
        af_models: mrparse_data.pfam.af_models,
        esm_models: mrparse_data.pfam.esm_models
    },
})
