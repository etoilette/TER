import spot
import buddy
import argparse

class my_twa:
    def __init__(self, aut, states_labels, outputs):
        self.twa = aut
        self.outputs = outputs
        self.states_labels = states_labels
        self.current_state = 0
    
    def step(self, inputs):
        #print(f" Moving with {inputs}")
        for edge in self.twa.out(self.current_state):
            #print(f"tries to go to {edge.dst}")
            if satisfies_bdd(inputs, edge.cond, self.outputs):
                self.current_state = edge.dst
                #print(f"Found an edge to {edge.dst}")
                return compute_outputs(self.twa, edge.cond, self.outputs)
        #print("found no edge")
        return compute_outputs(self.twa, self.states_labels[self.current_state], self.outputs)

def satisfies_bdd (valuation, bdd, outputs):
    if bdd == buddy.bddfalse:
        return False
    if bdd == buddy.bddtrue:
        return True
    top = buddy.bdd_var(bdd)
    low = buddy.bdd_low(bdd)
    high = buddy.bdd_high(bdd)
    if top in outputs:
        return satisfies_bdd(valuation, low, outputs) or satisfies_bdd(valuation, high, outputs)
    if valuation[top]:
        return satisfies_bdd(valuation, high, outputs)
    else:
        return satisfies_bdd(valuation, low, outputs)

def get_output_name(var, outputs):
    for name_out, var_out in outputs.items():
        if var_out == var:
            return name_out
    raise Exception(f"Variable {var_out} not found")

def compute_outputs(aut, bdd, outputs):
    resu = [False]*len(aut.ap())
    bdd_temp = bdd
    while bdd_temp != buddy.bddtrue:
        top = buddy.bdd_var(bdd_temp)
        low = buddy.bdd_low(bdd_temp)
        high = buddy.bdd_high(bdd_temp)
        if low == buddy.bddfalse:
            if top in outputs:
                resu[top] = True
            bdd_temp = high
            continue
        if high == buddy.bddfalse:
            bdd_temp = low
            continue
        raise Exception("Not a cube!")
    return resu

# Parsing of input

pn_text = """PetriNetType: 0
NumPlaces: 5
NumTransitions: 4
InitMark: [2,2,0,0,0]
<Transitions>
(0,"a"): [0],[2]|(i1;)
(1,"b"): [1],[3]|(;o2)
(2,"c"): [2,3],[4]|(i2;o1 & -o2)
(3,"x"): [0,1],[2,3]|(;)
</Transitions>
<Place-Output>
0:(o3)
1:(-o4 & o2)
</Place-Output>"""

def get_io_name(strg):
    if strg[0] == '!':
        return strg[1:]
    else:
        return strg

def parse_petri_net(pn_text):
    pn = dict()
    pn["Map"] = dict()
    pn["Transitions"] = dict()
    pn["Inputs"] = set()
    pn["Outputs"] = set()
    lines = pn_text.split("\n")

    i=0
    while lines[i] != "<Transitions>":
        i += 1
    i0 = i +1
    while lines[i] != "</Transitions>":
        i += 1
    init_marking = tuple([int(m) for m in ((lines[i0-2]).split(": ")[1][1:-1]).split(',')])
    pn["Initial"] = init_marking
    transitions_text = lines[i0:i]
    for transition_text in transitions_text:
        parts = transition_text.split(": ")
        name_texts = parts[0].split(",")
        id = int(name_texts[0][1:])
        name = name_texts[1][1:-2]
        pn["Map"][id] = name
        pn["Transitions"][name] = dict()
        defs = parts[1].split("|")
        defs_places = defs[0].split("],[")
        preplaces = [int(i) for i in defs_places[0][1:].split(",")]
        postplaces = [int(i) for i in defs_places[1][:-1].split(",")]
        pn["Transitions"][name]["pre"] = preplaces
        pn["Transitions"][name]["post"] = postplaces
        defs_io = defs[1][1:-1].split(";")
        inputs = defs_io[0]
        outputs = defs_io[1]
        inputs_names = [get_io_name(e) for e in inputs.split(" & ") if e != '']
        outputs_names = [get_io_name(e) for e in outputs.split(" & ") if e != '']
        pn["Transitions"][name]["inputs"] = inputs
        pn["Transitions"][name]["outputs"] = outputs
        pn["Inputs"].update(inputs_names)
        pn["Outputs"].update(outputs_names)
    
    pn["PlaceIO"] = dict()
    while lines[i] != "<Place-Output>":
        i += 1
    i0 = i +1
    while lines[i] != "</Place-Output>":
        i += 1
    if i > i0 +1:
        place_out_text = lines[i0:i]
        for placeio_text in place_out_text:
            placeio = placeio_text.split(':')
            outputs = placeio[1][1:-1]
            pn["PlaceIO"][int(placeio[0])] = outputs
            outputs_names = [get_io_name(e) for e in outputs.split(" & ") if e != '']
            pn["Outputs"].update(outputs_names)
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
    lines = hda_text.split("\n")
    hda = dict()
    i_start = 0
    while i_start<len(lines):
        id = lines[i_start].split(" : ")
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

        i_end = i_start+2
        while i_end < len(lines) and ":" not in lines[i_end]:
            i_end += 1
        transitions = lines[i_start+2:i_end-1]
        if transitions != []:
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
        i_start = i_end
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
            #print(f"{mset} is subsummed by {mset2}")
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

def get_literal(strg, ap):
    if strg[0] == "!":
        return buddy.bdd_not(ap[strg[1:]])
    else:
        return ap[strg]

def remove_space(strg):
    resu = ""
    for char in strg:
        if char !=' ':
            resu = resu + char
    return resu

def get_string_bdd(strg, ap):
    clauses = remove_space(strg).split('&')
    resu = buddy.bddtrue
    for clause in clauses:
        literals = clause.split('|')
        resu_clause = buddy.bddfalse
        for literal in literals:
            resu_clause = resu_clause | get_literal(literal, ap)
        resu = resu & resu_clause
    return resu

def get_marking_output(marking, pn, ap):
    resu = buddy.bddtrue
    #print(f"getting outputs for {marking} :")
    for state in range(len(marking)):
        if state in pn["PlaceIO"].keys() and marking[state] > 0:
            #print(f"getting bdd for state {state}: {pn["PlaceIO"][state]}")
            resu = resu & get_string_bdd(pn["PlaceIO"][state], ap)
    return resu

def get_state(marking, states, csg, pn, ap, states_label):
    if marking in states:
        state = states[marking]
    else:
        label = get_marking_output(marking, pn, ap)
        state = csg.new_state()
        states[marking] = state
        states_label[state] = label
    return state

class Inputs_Error (Exception):
    pass

class Outputs_Error (Exception):
    def __init__(self, message):            
        # Call the base class constructor with the parameters it needs
        super().__init__(message)

def get_io(marking_source, mset, pn, ap):
    resu = get_marking_output(marking_source, pn, ap)
    for t in mset.keys():
        inputs = pn["Transitions"][t]["inputs"]
        outputs = pn["Transitions"][t]["outputs"]
        #print(f"In mset: {mset}, inputs: \"{inputs}\" and outputs: \"{outputs}\"")
        if inputs:
            bdd_inputs = get_string_bdd(inputs, ap)
            resu = resu & bdd_inputs
            if resu == buddy.bddfalse:
                raise Inputs_Error
        if outputs:
            bdd_outputs = get_string_bdd(outputs, ap)
            resu = resu & bdd_outputs
            if resu == buddy.bddfalse:
                raise Outputs_Error(f"Incompatible outputs in {marking_source} when firing {mset}")
            resu = resu & get_string_bdd(outputs, ap)
    return resu

def add_transition(csg, pn, ap, states_label, edges_labels, marking_source, mset, states):
    marking_target = compute_marking_transi(marking_source, mset, pn)
    state_source = get_state(marking_source, states, csg, pn, ap, states_label)
    state_target = get_state(marking_target, states, csg, pn, ap, states_label)
    try:
        label = get_io(marking_source, mset, pn, ap)
    except Inputs_Error:
        return marking_target
    if not (marking_source, marking_target) in edges_labels.keys():
        edge = csg.new_edge(state_source, state_target, label)
        edges_labels[(marking_source, marking_target)] = [label]
        #print(f"added new edge from {marking_source} (state {state_source}) to {marking_target} (state {state_target}) by {mset}")
    elif not any(existing_label == label for existing_label in edges_labels[(marking_source, marking_target)]):
        edge = csg.new_edge(state_source, state_target, label)
        edges_labels[(marking_source, marking_target)].append(label)
        #print(f"added extra edge from {marking_source} (state {state_source}) to {marking_target} (state {state_target}) by {mset}")
    return marking_target

def add_many_transitions(csg, pn, ap, states_label, edges_labels, marking_source, concset, states, markings_avoidable):
    #print(f"-> computing transitions from {marking_source}x{concset}")
    # get_state(marking_source, states, csg, pn, ap)
    set_keys = set(concset.keys())
    for mset in sub_multi_set(concset, set_keys):
        if len(mset) != 0:
            marking_target = add_transition(csg, pn, ap, states_label, edges_labels, marking_source, mset, states)
            concset_comp = compute_comp_mset(mset, concset)
            #print(f"- finding {mset} to go from {marking_source} to {marking_target} with {concset_comp} remaining to do")
            #print(f"=== Avoiding consets {markings_avoidable}")
            if not marking_target in markings_avoidable:
                markings_avoidable.append(marking_target)
                #print(f"adding {mset} to avoidable concsets")
                #print(f"=== Now avoiding consets {markings_avoidable}")
                add_many_transitions(csg, pn, ap, states_label, edges_labels, marking_target, concset_comp, states, markings_avoidable)
                #print(f"-< finishing computing transitions from {marking_target}x{concset_comp}")
                #print(f"-> returning to computing transitions from {marking_source}x{concset}")
    return 1

""" This function takes an HDA in MAXCELL form that represents the petri net pn in inputs and produces its concurentstep graph associated"""
def hdatocsg (hda, pn) :
    bdict = spot.make_bdd_dict()
    csg = spot.make_twa_graph(bdict)

    ap = dict()
    for ap_in in pn["Inputs"]:
        ap_var = buddy.bdd_ithvar(csg.register_ap(ap_in))
        ap[ap_in] = ap_var
    
    outputs = set()
    for ap_out in pn["Outputs"]:
        ap_int = csg.register_ap(ap_out)
        ap_var = buddy.bdd_ithvar(ap_int)
        ap[ap_out] = ap_var
        outputs.add(ap_int)
        spot.set_synthesis_outputs(csg, ap_var)

    #states keeps track of the association marking <=> states of the automatons
    states = dict()

    edges_labels = dict()
    states_label = dict()

    initial_state = get_state(pn["Initial"], states, csg, pn, ap, states_label)
    csg.set_init_state(initial_state)
    for maxcell in hda.keys() :
        concset = hda[maxcell]["Concset"]
        marking = compute_marking_pre(hda[maxcell]["Marking"], concset , pn)
        #print(f"-- New MAXCELL {marking}x{concset}")
        markings_avoidable = compute_markings_avoidable(marking, hda[maxcell]["Transitions"], pn)
        add_many_transitions(csg, pn, ap, states_label, edges_labels, marking, concset, states, markings_avoidable)
    return my_twa(csg, states_label, outputs)

def dict_to_vect(dict1, dict2):
    resu = [False]*len(dict2)
    for ap in dict1:
        resu[dict2[ap]] = dict1[ap]
    return resu

if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    parser.add_argument("hda_file")
    parser.add_argument("pn_file")
    args = parser.parse_args()
    with open(args.hda_file) as hda_file, open(args.pn_file) as pn_file:
        pn = parse_petri_net(pn_file.read())
        #print(pn)
        hda = parse_maxcell_hda(hda_file.read(), pn)
        #print(hda)
        csg = hdatocsg(hda, pn)
        #print(csg.states_labels)
        print(csg.twa.to_str('hoa'))
        