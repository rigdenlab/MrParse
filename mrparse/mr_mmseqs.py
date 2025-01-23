# CODE ADAPTED FROM: OpenFold

import requests
import os
import time
import tarfile
import random
from typing import Tuple, List
import logging
logger = logging.getLogger(__name__)


def run(x, prefix, use_env=True, use_filter=True,
        use_templates=False, filter=None, use_pairing=False,
        host_url="https://a3m.mmseqs.com") -> Tuple[List[str], List[str]]:
  submission_endpoint = "ticket/pair" if use_pairing else "ticket/msa"

  def submit(seqs, mode, N=101):
    n, query = N, ""
    for seq in seqs:
      query += f">{n}\n{seq}\n"
      n += 1

    #print(f'requests.post({host_url}/{submission_endpoint})',{'q':query,'mode': mode})
    res = requests.post(f'{host_url}/{submission_endpoint}', data={'q':query,'mode': mode})
    try:
      out = res.json()
    except ValueError:
      logger.error(f"Server didn't reply with json: {res.text}")
      out = {"status":"ERROR"}
    return out

  def status(ID):
    res = requests.get(f'{host_url}/ticket/{ID}')
    try:
      out = res.json()
    except ValueError:
      logger.error(f"Server didn't reply with json: {res.text}")
      out = {"status":"ERROR"}
    return out

  def download(ID, path):
    res = requests.get(f'{host_url}/result/download/{ID}')
    with open(path,"wb") as out: out.write(res.content)

  # process input x
  seqs = [x] if isinstance(x, str) else x

  # compatibility to old option
  if filter is not None:
    use_filter = filter

  # setup mode
  if use_filter:
    mode = "env" if use_env else "all"
  else:
    mode = "env-nofilter" if use_env else "nofilter"

  if use_pairing:
    mode = ""
    use_templates = False
    use_env = False

  # define path
  #path = f"{prefix}_{mode}"
  path = prefix
  if not os.path.isdir(path): os.mkdir(path)

  # call mmseqs2 api
  tar_gz_file = f'{path}/out.tar.gz'
  N,REDO = 101,True

  # deduplicate and keep track of order
  seqs_unique = []
  #TODO this might be slow for large sets
  [seqs_unique.append(x) for x in seqs if x not in seqs_unique]
  Ms = [N + seqs_unique.index(seq) for seq in seqs]
  # lets do it!
  if not os.path.isfile(tar_gz_file):
    TIME_ESTIMATE = 150 * len(seqs_unique)
    while REDO:
        # Resubmit job until it goes through
        out = submit(seqs_unique, mode, N)
        while out["status"] in ["UNKNOWN", "RATELIMIT"]:
            sleep_time = 5 + random.randint(0, 5)
            logger.error(f"Sleeping for {sleep_time}s. Reason: {out['status']}")
            # resubmit
            time.sleep(sleep_time)
            out = submit(seqs_unique, mode, N)

        if out["status"] == "ERROR":
            raise Exception(f'MMseqs2 API is giving errors. Please confirm your input is a valid protein sequence. If error persists, please try again an hour later.')

        if out["status"] == "MAINTENANCE":
            raise Exception(f'MMseqs2 API is undergoing maintenance. Please try again in a few minutes.')

        # wait for job to finish
        ID,TIME = out["id"],0
        while out["status"] in ["UNKNOWN","RUNNING","PENDING"]:
            t = 5 + random.randint(0,5)
            logger.error(f"Sleeping for {t}s. Reason: {out['status']}")
            time.sleep(t)
            out = status(ID)
            if out["status"] == "RUNNING":
                TIME += t

        if out["status"] == "COMPLETE":
            REDO = False

        if out["status"] == "ERROR":
            REDO = False
            raise Exception(f'MMseqs2 API is giving errors. Please confirm your input is a valid protein sequence. If error persists, please try again an hour later.')

        # Download results
        download(ID, tar_gz_file)

  # prep list of a3m files
  if use_pairing:
    a3m_files = [f"{path}/pair.a3m"]
  else:
    a3m_files = [f"{path}/uniref.a3m"]
    if use_env: a3m_files.append(f"{path}/bfd.mgnify30.metaeuk30.smag30.a3m")

  # extract a3m files
  if any(not os.path.isfile(a3m_file) for a3m_file in a3m_files):
    with tarfile.open(tar_gz_file) as tar_gz:
      tar_gz.extractall(path)

  # templates
  if use_templates:
    templates = {}
    #print("seq\tpdb\tcid\tevalue")
    for line in open(f"{path}/pdb70.m8","r"):
      p = line.rstrip().split()
      M,pdb,qid,e_value = p[0],p[1],p[2],p[10]
      M = int(M)
      if M not in templates: templates[M] = []
      templates[M].append(pdb)
      #if len(templates[M]) <= 20:
      #  print(f"{int(M)-N}\t{pdb}\t{qid}\t{e_value}")

    template_paths = {}
    for k,TMPL in templates.items():
      TMPL_PATH = f"{prefix}_{mode}/templates_{k}"
      if not os.path.isdir(TMPL_PATH):
        os.mkdir(TMPL_PATH)
        TMPL_LINE = ",".join(TMPL[:20])
        os.system(f"curl -s https://a3m-templates.mmseqs.com/template/{TMPL_LINE} | tar xzf - -C {TMPL_PATH}/")
        os.system(f"cp {TMPL_PATH}/pdb70_a3m.ffindex {TMPL_PATH}/pdb70_cs219.ffindex")
        os.system(f"touch {TMPL_PATH}/pdb70_cs219.ffdata")
      template_paths[k] = TMPL_PATH

  # gather a3m lines
  a3m_lines = {}
  for a3m_file in a3m_files:
    update_M,M = True,None
    for line in open(a3m_file,"r"):
      if len(line) > 0:
        if "\x00" in line:
          line = line.replace("\x00","")
          update_M = True
        if line.startswith(">") and update_M:
          M = int(line[1:].rstrip())
          update_M = False
          if M not in a3m_lines: a3m_lines[M] = []
        a3m_lines[M].append(line)

  # return results

  a3m_lines = ["".join(a3m_lines[n]) for n in Ms]

  if use_templates:
    template_paths_ = []
    for n in Ms:
      if n not in template_paths:
        template_paths_.append(None)
        #print(f"{n-N}\tno_templates_found")
      else:
        template_paths_.append(template_paths[n])
    template_paths = template_paths_


  return (a3m_lines, template_paths) if use_templates else a3m_lines


if __name__ == "__main__":
  import argparse
  parser = argparse.ArgumentParser(description='MMseqs2')
  parser.add_argument('--prefix', default="tmp_mmseqs2", help='prefix for output files')
  parser.add_argument('--use_env', default=True, type=bool, help='use env or not')
  parser.add_argument('--use_filter', default=True, type=bool, help='use filter or not')
  parser.add_argument('--use_templates', default=False, type=bool, help='use templates or not')
  parser.add_argument('--use_pairing', default=False, type=bool, help='use pairing or not')
  parser.add_argument('--host_url', default="https://a3m.mmseqs.com", help='host url')
  parser.add_argument('seqs', type=str, nargs='+', help='protein sequences')
  args = parser.parse_args()

  for seq in args.seqs:
    a3m_lines, template_paths = run(seq, args.prefix, args.use_env, args.use_filter, args.use_templates, use_pairing=args.use_pairing, host_url=args.host_url)
    print(a3m_lines)
    print(template_paths)
    print("")