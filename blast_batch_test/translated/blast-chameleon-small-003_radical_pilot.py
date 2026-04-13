#!/usr/bin/env python
import radical.pilot as rp

session = rp.Session()
try:
    pmgr = rp.PilotManager(session=session)
    tmgr = rp.TaskManager(session=session)

    pilot = pmgr.submit_pilots(rp.PilotDescription({
        'resource'     : 'local.localhost',
        'cores'        : 4,
        'runtime'      : 30,
        'exit_on_error': False
    }))
    tmgr.add_pilots(pilot)
    pilot.wait(rp.PMGR_ACTIVE)

    # ---- Level 0: split_fasta_ID000001 (no parents) ----
    td_split_fasta_ID000001 = rp.TaskDescription({
        'uid'           : 'split_fasta_ID000001',
        'executable'    : 'split_fasta',
        'arguments'     : ['./split_fasta', '5', 'small.fasta'],
        'cores_per_rank': 1,
        'input_staging' : [],
        'output_staging': [
            {'source': 'task:///small.fasta.0',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.0',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.1',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.1',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.2',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.2',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.3',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.3',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.4',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.4',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.5',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.5',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.6',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.6',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.7',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.7',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.8',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.8',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.9',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.9',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.10',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.10',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.11',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.11',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.12',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.12',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.13',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.13',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.14',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.14',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.15',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.15',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.16',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.16',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.17',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.17',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.18',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.18',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.19',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.19',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.20',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.20',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.21',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.21',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.22',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.22',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.23',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.23',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.24',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.24',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.25',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.25',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.26',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.26',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.27',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.27',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.28',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.28',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.29',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.29',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.30',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.30',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.31',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.31',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.32',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.32',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.33',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.33',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.34',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.34',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.35',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.35',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.36',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.36',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.37',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.37',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.38',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.38',
             'action': rp.COPY},
            {'source': 'task:///small.fasta.39',
             'target': 'pilot:///split_fasta_ID000001_small.fasta.39',
             'action': rp.COPY}
        ]
    })
    tasks_0 = tmgr.submit_tasks([td_split_fasta_ID000001])
    tmgr.wait_tasks([t.uid for t in tasks_0])

    # ---- Level 1: blastall_ID000002 to blastall_ID000041 ----
    blast_tasks = []
    for i in range(42):
        task_id = f'blastall_ID00000{i+2}'
        blast_task = rp.TaskDescription({
            'uid'           : task_id,
            'executable'    : 'blastall',
            'arguments'     : [f'./blastall', '-p', 'blastn', '-d', 'nt/nt', '-i', f'small.fasta.{i}', '-o', f'small.fasta.{i}.out', '2>', f'small.fasta.{i}.err'],
            'cores_per_rank': 1,
            'input_staging' : [
                {'source': 'pilot:///split_fasta_ID000001_small.fasta.' + str(i),
                 'target': 'task:///small.fasta.' + str(i),
                 'action': rp.LINK},
                {'source': 'pilot:///blastall',
                 'target': 'task:///blastall',
                 'action': rp.LINK},
                {'source': 'client:///data/nt',
                 'target': 'task:///nt',
                 'action': rp.TRANSFER}
            ],
            'output_staging': [
                {'source': 'task:///small.fasta.' + str(i) + '.out',
                 'target': 'pilot:///blastall_ID00000' + str(i+2) + '_small.fasta.' + str(i) + '.out',
                 'action': rp.COPY},
                {'source': 'task:///small.fasta.' + str(i) + '.err',
                 'target': 'pilot:///blastall_ID00000' + str(i+2) + '_small.fasta.' + str(i) + '.err',
                 'action': rp.COPY}
            ]
        })
        blast_tasks.append(blast_task)
    tasks_1 = tmgr.submit_tasks(blast_tasks)
    tmgr.wait_tasks([t.uid for t in tasks_1])

    # ---- Level 2: cat_blast_ID000042, cat_ID000043 ----
    td_cat_blast_ID000042 = rp.TaskDescription({
        'uid'           : 'cat_blast_ID000042',
        'executable'    : 'cat_blast',
        'arguments'     : ['./cat_blast', 'None'] + [f'small.fasta.{i}.out' for i in range(40)],
        'cores_per_rank': 1,
        'input_staging' : [],
        'output_staging': [
            {'source': 'task:///None',
             'target': 'client:///data/None',
             'action': rp.TRANSFER}
        ]
    })
    td_cat_ID000043 = rp.TaskDescription({
        'uid'           : 'cat_ID000043',
        'executable'    : 'cat',
        'arguments'     : ['cat'] + [f'small.fasta.{i}.err' for i in range(40)] + ['>', 'None.err'],
        'cores_per_rank': 1,
        'input_staging' : [],
        'output_staging': [
            {'source': 'task:///None.err',
             'target': 'client:///data/None.err',
             'action': rp.TRANSFER}
        ]
    })
    tasks_2 = tmgr.submit_tasks([td_cat_blast_ID000042, td_cat_ID000043])
    tmgr.wait_tasks([t.uid for t in tasks_2])

finally:
    session.close()