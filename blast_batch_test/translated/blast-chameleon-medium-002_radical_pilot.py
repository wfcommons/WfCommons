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
        'executable'    : './split_fasta',
        'arguments'     : ['100', 'medium.fasta'],
        'cores_per_rank': 1,
        'output_staging': [
            {'source': 'task:///medium.fasta.0',
             'target': 'pilot:///medium.fasta.0',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.1',
             'target': 'pilot:///medium.fasta.1',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.2',
             'target': 'pilot:///medium.fasta.2',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.3',
             'target': 'pilot:///medium.fasta.3',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.4',
             'target': 'pilot:///medium.fasta.4',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.5',
             'target': 'pilot:///medium.fasta.5',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.6',
             'target': 'pilot:///medium.fasta.6',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.7',
             'target': 'pilot:///medium.fasta.7',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.8',
             'target': 'pilot:///medium.fasta.8',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.9',
             'target': 'pilot:///medium.fasta.9',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.10',
             'target': 'pilot:///medium.fasta.10',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.11',
             'target': 'pilot:///medium.fasta.11',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.12',
             'target': 'pilot:///medium.fasta.12',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.13',
             'target': 'pilot:///medium.fasta.13',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.14',
             'target': 'pilot:///medium.fasta.14',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.15',
             'target': 'pilot:///medium.fasta.15',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.16',
             'target': 'pilot:///medium.fasta.16',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.17',
             'target': 'pilot:///medium.fasta.17',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.18',
             'target': 'pilot:///medium.fasta.18',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.19',
             'target': 'pilot:///medium.fasta.19',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.20',
             'target': 'pilot:///medium.fasta.20',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.21',
             'target': 'pilot:///medium.fasta.21',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.22',
             'target': 'pilot:///medium.fasta.22',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.23',
             'target': 'pilot:///medium.fasta.23',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.24',
             'target': 'pilot:///medium.fasta.24',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.25',
             'target': 'pilot:///medium.fasta.25',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.26',
             'target': 'pilot:///medium.fasta.26',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.27',
             'target': 'pilot:///medium.fasta.27',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.28',
             'target': 'pilot:///medium.fasta.28',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.29',
             'target': 'pilot:///medium.fasta.29',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.30',
             'target': 'pilot:///medium.fasta.30',
             'action': rp.COPY},
            {'source': 'task:///medium.fasta.31',
             'target': 'pilot:///medium.fasta.31',
             'action': rp.COPY}
        ]
    })
    tasks_0 = tmgr.submit_tasks([td_split_fasta_ID000001])
    tmgr.wait_tasks([t.uid for t in tasks_0])

    # ---- Level 1: blastall_ID000002 to blastall_ID000301 ----
    blast_tasks = []
    for i in range(302):
        td_blastall = rp.TaskDescription({
            'uid'           : f'blastall_ID000{i:03}',
            'executable'    : './blastall',
            'arguments'     : ['-p', 'blastn', '-d', 'nt/nt', '-i', f'medium.fasta.{i}', '-o', f'medium.fasta.{i}.out', '2>', f'medium.fasta.{i}.err'],
            'cores_per_rank': 1,
            'input_staging' : [
                {'source': 'pilot:///split_fasta_ID000001_medium.fasta',
                 'target': f'task:///medium.fasta.{i}',
                 'action': rp.LINK}
            ],
            'output_staging': [
                {'source': f'task:///medium.fasta.{i}.out',
                 'target': f'pilot:///medium.fasta.{i}.out',
                 'action': rp.COPY},
                {'source': f'task:///medium.fasta.{i}.err',
                 'target': f'pilot:///medium.fasta.{i}.err',
                 'action': rp.COPY}
            ]
        })
        blast_tasks.append(td_blastall)

    tasks_1 = tmgr.submit_tasks(blast_tasks)
    tmgr.wait_tasks([t.uid for t in tasks_1])

    # ---- Level 2: cat_blast_ID000302 ----
    td_cat_blast_ID000302 = rp.TaskDescription({
        'uid'           : 'cat_blast_ID000302',
        'executable'    : './cat_blast',
        'arguments'     : ['None'] + [f'medium.fasta.{i}.out' for i in range(302)],
        'cores_per_rank': 1,
        'input_staging' : [
            {'source': 'pilot:///split_fasta_ID000001_medium.fasta',
             'target': 'task:///cat_blast',
             'action': rp.LINK}
        ],
        'output_staging': [
            {'source': 'task:///None',
             'target': 'client:///data/cat_blast.out',
             'action': rp.TRANSFER}
        ]
    })
    tasks_2 = tmgr.submit_tasks([td_cat_blast_ID000302])
    tmgr.wait_tasks([t.uid for t in tasks_2])

    # ---- Level 3: cat_ID000303 ----
    td_cat_ID000303 = rp.TaskDescription({
        'uid'           : 'cat_ID000303',
        'executable'    : 'cat',
        'arguments'     : ['cat'] + [f'medium.fasta.{i}.err' for i in range(302)],
        'cores_per_rank': 1,
        'output_staging': [
            {'source': 'task:///None.err',
             'target': 'client:///data/cat.err',
             'action': rp.TRANSFER}
        ]
    })
    tasks_3 = tmgr.submit_tasks([td_cat_ID000303])
    tmgr.wait_tasks([t.uid for t in tasks_3])

finally:
    session.close()