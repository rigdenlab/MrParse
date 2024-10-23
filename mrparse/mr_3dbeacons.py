# Putting all redundant 3Dbeacons code here incase it improves and we want to reimplement it

def run_3dbeacons_alphafold_api(seq_info, max_hits=10):
    url = 'https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/sequence/search'
    headers = {
    'accept': 'application/json',
    'Content-Type': 'application/json'
    }
    data = {
    'sequence': seq_info.sequence,
    }

    response = requests.post(url, headers=headers, json=data)
    job_id = json.loads(response.text)['job_id']

    url = 'https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/sequence/result'
    params = {
        'job_id': job_id,
    }
    headers = {
        'accept': 'application/json'
    }

    response = requests.get(url, params=params, headers=headers)
    data = json.loads(response.text)

    while 'message' in data and data['message'] == 'Search in progress, please try after sometime!':
        time.sleep(5)
        url = 'https://www.ebi.ac.uk/pdbe/pdbe-kb/3dbeacons/api/sequence/result'
        params = {
            'job_id': job_id,
        }
        headers = {
            'accept': 'application/json'
        }

        response = requests.get(url, params=params, headers=headers)
        data = json.loads(response.text)
    
    if 'message' in data and data['message'] == "Error in submitting the job, please retry!":
        logger.error("3D-Beacons API error: %s" % data['message'])
        return []

    alphafold_structures = []
    for entry in data:
        hsp_data = entry['hit_hsps']
        for structure in entry['summary']['structures']:
            provider = structure['summary']['provider']
            if provider == 'AlphaFold DB':
                alphafold_structures.append([hsp_data, structure])

    return alphafold_structures


def _find_api_hits(json_list, input_sequence, max_hits=10):
    hitDict = OrderedDict()
    for structure in json_list:
        hsp_data = structure[0][0]
        data = structure[1]['summary']

        model_id = data['entities'][0]['identifier']
        base_url = "http://www.uniprot.org/uniprot/"
        current_url = base_url + model_id + ".fasta"
        response = requests.post(current_url)
        model_sequence = ''.join(response.text)

        with open("downloaded_sequences.fasta", "w") as f:
            f.write(f"{model_id}\n{model_sequence}\n")

    from mrparse.searchDB import phmmer
    run_phmmer(input_sequence, pdb_seqdb="downloaded_sequences.fasta", dblvl="model")
    phr = phmmer()
    logfile = "phmmer_model.log"
    plog = open(logfile, "r")
    phmmerALNLog = plog.readlines()
    plog.close()
    phr.logfile = logfile
    phr.getPhmmerAlignments(targetSequence=input_sequence, phmmerALNLog=phmmerALNLog, PDBLOCAL=None, DB="simple", seqMetaDB=None)

    for hitname in phr.resultsDict.keys():
        try:
            sh = SequenceHit()
            sh.rank = phr.resultsDict[hitname].rank

            sh.pdb_id = phr.resultsDict[hitname].afdbName
            sh.evalue = phr.resultsDict[hitname].evalue
            sh.query_start = phr.resultsDict[hitname].tarRange[0]
            sh.query_stop = phr.resultsDict[hitname].tarRange[1]
            sh.hit_start = int(phr.resultsDict[hitname].alnRange[0])
            sh.hit_stop = int(phr.resultsDict[hitname].alnRange[1])
            sh.target_alignment = phr.resultsDict[hitname].targetAlignment
            sh.alignment = phr.resultsDict[hitname].alignment

            hstart = sh.hit_start
            hstop = sh.hit_stop
            qstart, qstop = sh.query_start, sh.query_stop

            sh.query_start = qstart
            sh.query_stop = qstop
            sh.hit_start = hstart
            sh.hit_stop = hstop

            seq_ali = zip(range(qstart, qstop), sh.alignment)
            sh.seq_ali = [x[0] for x in seq_ali if x[1] != '-']

            local, overall = phr.resultsDict[hitname].localSEQID, phr.resultsDict[hitname].overallSEQID
            sh.local_sequence_identity = np.round(local)
            sh.overall_sequence_identity = np.round(overall)

            sh.search_engine = "3D-Beacons"
            sh.model_url = data['model_url']
            sh.data_created = data['created']

            if sh.rank <= max_hits:
                    hitDict[hitname] = sh

        except Exception:
            logger.debug(f"Issue with target {data['entities'][0]['identifier']}")
    os.unlink("downloaded_sequences.fasta")
    return hitDict