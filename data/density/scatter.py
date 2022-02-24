import numpy as np


klms = {
    "den-sym": [0.39269908169, 0, -0.09766608, 0, 0, 0, 0, 0, 0, 0],
    "den-asym": [0.39269908169, 0.05200629, -0.2021978, 0, 0, 0, 0, 0, 0, 0],
    "den-tet": [],
    "den-db": [],
    "den-high": [0, 0, 0, 0.004444906415833378, 0.003477837670425714, 0.005378876304115907, 0.005783583661113503, -0.002305670990235569, -0.004974910518199871, 0.003351084537350686],
    "den-in": [],
    "den-out": [],
}

def get_text(params):
    return """0, 3
120
5
1000
6000
0.00006464182, 0.00012928364, -0.00012928364
1.0
{}, {}, {}, {}, {}, {}, {}, {}, {}, {}
0.78539816339, 0.125, 0, 1, 1, 1, 1, 1, 1, 1
-0.78539816339, -0.125, -0.25, -1, -1, -1, -1, -1, -1, -1
1e-2, 1e-5""".format(params[0], params[1], params[2], params[3], params[4],
    params[5], params[6], params[7], params[8], params[9])

for name, params in klms.items():
    with open("../../staged/{}.txt".format(name), 'w') as f:
        f.write(get_text(params))
