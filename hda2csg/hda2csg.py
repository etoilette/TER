import spot
import buddy
import argparse

# Parsing of input

pn_text = """PetriNetType: 0
NumPlaces: 5
NumTransitions: 4
InitMark: [2,2,0,0,0]
<Transitions>
(0,"a"): [0],[2]
(1,"b"): [1],[3]
(2,"c"): [2,3],[4]
(3,"x"): [0,1],[2,3]
</Transitions>"""

def parse_petri_net(pn_text):
    lines = pn_text.split("\n")
    i=0
    while lines[i] != "<Transitions>":
        i += 1
    i0 = i +1
    while lines[i] != "</Transitions>":
        i += 1
    transitions_text = lines[i0:i]
    pn = dict()
    pn["Map"] = dict()
    pn["Transitions"] = dict()
    for transition_text in transitions_text:
        parts = transition_text.split(": ")
        name_texts = parts[0].split(",")
        id = int(name_texts[0][1:])
        name = name_texts[1][1:-2]
        pn["Map"][id] = name
        pn["Transitions"][name] = dict()
        defs = parts[1].split("],[")
        preplaces = [int(i) for i in defs[0][1:].split(",")]
        postplaces = [int(i) for i in defs[1][:-1].split(",")]
        pn["Transitions"][name]["pre"] = preplaces
        pn["Transitions"][name]["post"] = postplaces
    return pn

hda_text = """4 : [0,0,0,0,0] x [a,a,b,b]
edges:
1;()|(0,0,1,1);()|(2,2)|();10
1;()|(0,1);(0,1)|(2)|();17
1;(0,1)|();(0,1)|(3)|();46
0;()|(0,1);(0,1)|()|(3);46
1;(0,1)|(0,1);()|(2,3)|();49
0;()|(0,0,1,1);()|(2)|(3);49
1;(0,0,1,1)|();()|(3,3)|();68
0;(0,1)|(0,1);()|(3)|(3);68
0;()|(0,0,1,1);()|()|(3,3);68

17 : [0,0,0,0,0] x [a,b,c]
edges:
0;(2)|();(0,1)|()|(3);46
1;(0,1)|();(2)|(3)|();49
0;()|(0,1);(2)|()|(3);49
0;(0,1,2)|();()|(3)|(3);68
0;(2)|(0,1);()|()|(3,3);68

46 : [0,0,0,0,0] x [a,b,x]
edges:
1;(0,1)|(3);()|(2,3)|();49
1;()|(0,1);(3)|(2)|();49
1;(0,1)|();(3)|(3)|();68
0;()|(0,1);(3)|()|(3);68

10 : [0,0,0,0,0] x [c,c]
edges:
0;(2)|();(2)|()|(0,1);17
0;(2,2)|();()|()|(0,1,3);46
0;(2)|();(2)|()|(3);49
0;(2,2)|();()|()|(3,3);68

49 : [0,0,0,0,0] x [c,x]
edges:
0;(2)|();(3)|()|(3);68

68 : [0,0,0,0,0] x [x,x]
edges:
"""

def parse_maxcell_hda(hda_text, pn):
    Maxcells = hda_text.split("\n\n")
    hda = dict()
    for maxcell_text in Maxcells:
        lines = maxcell_text.split("\n")
        id = lines[0].split(" : ")
        name = int(id[0])
        hda[name] = dict()
        markconc = id[1].split(" x ")
        marking = tuple([int(i) for i in markconc[0][1:-1].split(",")])
        concset = dict()
        for a in markconc[1][1:-1].split(","):
            if a in concset.keys():
                concset[a] += 1
            else:
                concset[a] = 1
        hda[name]["Marking"] = marking
        hda[name]["Concset"] = concset
        hda[name]["Transitions"] = dict()
        transitions = lines[2:]
        if transitions[0] != '':
            for transition_text in transitions:
                transition_parts = transition_text.split(";")
                maxcell_target = int(transition_parts[3])
                if maxcell_target not in hda[name]["Transitions"].keys():
                    hda[name]["Transitions"][maxcell_target] = []
                transition = dict()
                transition["Unstart"] = dict()
                transition["Terminate"] = dict()
                transition["Continued"] = dict()
                transition["Start"] = dict()
                unst_term = transition_parts[1].split("|")
                unstart_list = unst_term[0][1:-1].split(",")
                continued_list = unst_term[0][1:-1].split(",")
                terminate_list = unst_term[1][1:-1].split(",")
                start_list = transition_parts[2].split("|")[1][1:-1].split(",")
                for a in unstart_list:
                    if a != '':
                        char = pn["Map"][int(a)]
                        if char in transition["Unstart"].keys():
                            transition["Unstart"][char] += 1
                        else:
                            transition["Unstart"][char] = 1
                for a in terminate_list:
                    if a != '':
                        char = pn["Map"][int(a)]
                        if char in transition["Terminate"].keys():
                            transition["Terminate"][char] += 1
                        else:
                            transition["Terminate"][char] = 1
                for a in continued_list:
                    if a != '':
                        char = pn["Map"][int(a)]
                        if char in transition["Continued"].keys():
                            transition["Continued"][char] += 1
                        else:
                            transition["Continued"][char] = 1
                for a in start_list:
                    if a != '':
                        char = pn["Map"][int(a)]
                        if char in transition["Start"].keys():
                            transition["Start"][char] += 1
                        else:
                            transition["Start"][char] = 1
                hda[name]["Transitions"][maxcell_target].append(transition)
    return hda

# Computing concurrent step graph

def compute_comp_mset(mset, bigmset):
    resu = dict()
    for key in bigmset.keys():
        if key in mset.keys():
            if mset[key] < bigmset[key]:
                resu[key] = bigmset[key] - mset[key]
        else:
            resu[key] = bigmset[key]
    return resu

def diff_mset(mset1, mset2):
    resu = dict()
    for key in mset1.keys():
        if key in mset2.keys():
            value = mset1[key] - mset2[key]
            if value > 0:
                resu[key] = value
        else:
            resu[key] = mset1[key]
    return resu

def leq_mset(mset1, mset2):
    for key in mset1.keys():
        if key not in mset2.keys() or mset1[key] > mset2[key]:
            return False
    return True

def is_subsumed_msets(mset, msets):
    for mset2 in msets:
        if leq_mset(mset2, mset):
            print(f"{mset} is subsummed by {mset2}")
            return True
    return False

def sub_multi_set(mset, keys):
    #print(keys)
    if len(keys) == 0:
        yield dict()
    else:
        key = keys.pop()
        #print(f"working on {key} for {mset[key]}")
        for sub_mset in sub_multi_set(mset, keys):
            for i in range(mset[key]+1):
                if i != 0:
                    sub_mset[key] = i
                elif key in sub_mset.keys():
                    del sub_mset[key]
                yield sub_mset.copy()
        keys.add(key)

def compute_marking_pre(marking, concset, pn) :
    new_marking = ()
    for place in range(len(marking)):
        mplace = marking[place]
        for transition in concset.keys():
            if place in pn["Transitions"][transition]["pre"] :
                mplace += concset[transition]
        new_marking = new_marking + (mplace,)
    # print(f"going from {marking} to {new_marking}")
    return new_marking

def compute_marking_transi(marking, concset, pn) :
    new_marking = ()
    for place in range(len(marking)):
        mplace = marking[place]
        for transition in concset.keys():
            if place in pn["Transitions"][transition]["pre"] :
                mplace -= concset[transition]
                if mplace < 0:
                    raise Exception("transition with not enough tokens")
            if place in pn["Transitions"][transition]["post"] :
                mplace += concset[transition]
        new_marking = new_marking + (mplace,)
    # print(f"going from {marking} to {new_marking} with {concset}")
    return new_marking

def compute_marking_maxcell_transi(marking, transitions, pn) :
    new_marking = ()
    for place in range(len(marking)):
        mplace = marking[place]
        for unstart in transitions["Unstart"].keys():
            if place in pn["Transitions"][unstart]["pre"] :
                mplace += transitions["Unstart"][unstart]
        for terminate in transitions["Terminate"].keys():
            if place in pn["Transitions"][terminate]["post"] :
                mplace += transitions["Terminate"][terminate]
        for start in transitions["Start"].keys():
            if place in pn["Transitions"][start]["pre"] :
                mplace -= transitions["Start"][start]
                if mplace < 0:
                    raise Exception("transition with not enough tokens")
        new_marking = new_marking + (mplace,)
    # print(f"going from {marking} to {new_marking} with {concset}")
    return new_marking

def has_face_covered(concset, transitions_list):
    set_keys = set(concset.keys())
    for subconcset in sub_multi_set(concset, set_keys):
        if not any(subconcset == transition["Continued"] for transition in transitions_list):
            return False
    return True

def compute_markings_avoidable(marking, transitions_dict, pn):
    resu = []
    for maxcell in transitions_dict:
        for transitions in transitions_dict[maxcell]:
            if has_face_covered(transitions["Continued"], transitions_dict[maxcell]):
                resu.append(compute_marking_maxcell_transi(marking, transitions, pn))
    return resu

def add_transition(csg, pn, marking_source, mset, states, p1):
    marking_target = compute_marking_transi(marking_source, mset, pn)
    """if mark ing_source in states:
        state_source = states[marking_source]
    else:
        state_source = csg.new_state()
        states[marking_source] = state_source
        print(f"new state {state_source} for source {marking_source}") """
    state_source = states[marking_source]
    if marking_target in states:
        state_target = states[marking_target]
    else:
        state_target = csg.new_state()
        states[marking_target] = state_target
        print(f"new state {state_target} for {marking_target}")
    if not any(e.dst == state_target for e in csg.out(state_source)):
        csg.new_edge(state_source, state_target, p1)
    print(f"added edge from {marking_source} (state {state_source}) to {marking_target} (state {state_target}) by {mset}")
    return marking_target

def add_many_transitions(csg, pn, marking_source, concset, states, markings_avoidable, p1):
    print(f"-> computing transitions from {marking_source}x{concset}")
    if marking_source in states:
        state = states[marking_source]
    else:
        state = csg.new_state()
        states[marking_source] = state
        print(f"new state {state} for source {marking_source}")
    set_keys = set(concset.keys())
    for mset in sub_multi_set(concset, set_keys):
        if len(mset) != 0:
            marking_target = add_transition(csg, pn, marking_source, mset, states, p1)
            concset_comp = compute_comp_mset(mset, concset)
            print(f"- finding {mset} to go from {marking_source} to {marking_target} with {concset_comp} remaining to do")
            print(f"=== Avoiding consets {markings_avoidable}")
            if not marking_target in markings_avoidable:
                markings_avoidable.append(marking_target)
                print(f"adding {mset} to avoidable concsets")
                print(f"=== Now avoiding consets {markings_avoidable}")
                add_many_transitions(csg, pn, marking_target, concset_comp, states, markings_avoidable, p1)
                print(f"-< finishing computing transitions from {marking_target}x{concset_comp}")
                print(f"-> returning to computing transitions from {marking_source}x{concset}")
    return 1

""" This function takes an HDA in MAXCELL form that represents the petri net pn in inputs and produces its concurentstep graph associated"""
def hdatocsg (hda, pn) :
    bdict = spot.make_bdd_dict()
    csg = spot.make_twa_graph(bdict)
    p1 = buddy.bdd_ithvar(csg.register_ap("p1"))

    #states keeps track of the association marking <=> states of the automatons
    states = dict()

    for maxcell in hda.keys() :
        concset = hda[maxcell]["Concset"]
        marking = compute_marking_pre(hda[maxcell]["Marking"], concset , pn)
        print(f"-- New MAXCELL {marking}x{concset}")
        markings_avoidable = compute_markings_avoidable(marking, hda[maxcell]["Transitions"], pn)
        add_many_transitions(csg, pn, marking, concset, states, markings_avoidable, p1)
    return csg


if __name__ == "__main__" :
    # Utiliser ARGPARSE pour les argumments !
    parser = argparse.ArgumentParser()
    parser.add_argument("hda_file")
    parser.add_argument("pn_file")
    args = parser.parse_args()
    with open(args.hda_file) as hda_file, open(args.pn_file) as pn_file:
        print("Hello\n")
        pn = parse_petri_net(pn_file.read())
        print(pn)
        hda = parse_maxcell_hda(hda_file.read(), pn)
        print(hda)
        csg = hdatocsg(hda, pn)
        print(csg.to_str('hoa'))
    # mset = {"a" : 2, "b" : 2}
    # keys = set(mset.keys())
    # for submset in sub_multi_set(mset, keys):
    #     print(submset)