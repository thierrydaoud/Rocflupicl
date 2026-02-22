#the props dict will have the property name as the key
#and an array with 3 elements, 0 is the datatype, 1 is bool for being a solution, 2 is bool for ghost prop
props = {}
arrayLens = {}
FinalStructs = None
Mapping = {}
ndim = -1


# five types of arrays: pos, real_sol, real_sol_wpos(real_sol + pos), int_sol, real_prop, int_prop, proj (projection). Each one defines one or more arrays in one or more structures.
# this allows us to swap layouts if we want to optimize memory accesses
# the third element of each component array specification, is a list of ways that the array can be refered to, ie y or ydot
layout = [
    {
        "StructName": "PPICLF_t_particlepos",
        "ppiclfArrayName" : "ppiclf_partpos",
        "ppiclfGhostArrayName" : "ppiclf_partpos_ghost",
        "arrays":
            [
                ("pos",  "pos", ["y"])
                ]
    },
    {
        "StructName": "PPICLF_t_particle",
        "ppiclfArrayName" : "ppiclf_parts",
        "ppiclfGhostArrayName" : "ppiclf_parts_ghost",
        "arrays":
            [
                ("y_real",  "real_sol", ["y"]),
                ("y_int",   "int_sol", ["y"]),
                ("ydot_real",   "real_sol_wpos", ["ydot"]),
                ("ydot_int",    "int_sol", ["ydot"]),
                ("ydotc_real", "real_sol_wpos", ["ydotc"]),
                ("ydotc_int", "int_sol", ["ydotc"]),
                ("rprop",   "real_prop", ["prop", "rprop"]),
                ("iprop",   "int_prop", ["prop", "iprop"]),
                ("feedback", "proj", ["feedbk", "proj", "feedback"]),
                ("y1_real", "real_sol_wpos", ["y1"]),
                ("y1_int", "int_sol", ["y1"]),
                ]
    }
]

INTERNAL_evaldtype = lambda args: "int" if "int" in args else "real"
INTERNAL_evalsolution = lambda args: "solution" in args
INTERNAL_evalghost = lambda args: "ghost" in args
INTERNAL_evaltemp = lambda args: "temp" in args
INTERNAL_evalproj = lambda args: "proj" in args

INTERNAL_evalargs = lambda args : [INTERNAL_evaldtype(args), INTERNAL_evalsolution(args), INTERNAL_evalghost(args), INTERNAL_evaltemp(args), INTERNAL_evalproj(args)]

def setvar(name, value):
    globals()[name] = value


def AddProp(name, *args):
    name = name.lower()
    if name in props.keys():
        raise RuntimeError("Property " + name + " redefined")
    props.update({name : INTERNAL_evalargs(args)})


# This is the function that actually decides what indexes to assign to each property in the array
# we could dump the results of this to a file, and have everything that uses fypp depend on it to ensure that we 
# never end up with a messed up mapping, but that shouldn't be needed
def GenerateMapping(**kwargs):
    global Mapping, arrayLens
    Mapping = {}
    for arrayName, arrayProps in kwargs.items():
        # print(arrayName)
        arrayLens[arrayName] = len(arrayProps)
        if len(arrayProps) == 0:
            # print("Skipping " + arrayName)
            continue
        # this sort puts all ghost properties at the begining of the array
        finalOrder = sorted(arrayProps, key=lambda x: 0 if x[1] else 1)
        ghostMapping = list(map(lambda x: x[1], finalOrder))
        if True in ghostMapping:
            lastGhostIndex = len(finalOrder) - ghostMapping[::-1].index(True)
        else:
            lastGhostIndex = 0
        # print(finalOrder)
        # print(lastGhostIndex)
        arrayMapping = { prop[0] : i + 1 for i, prop in enumerate(finalOrder)}
        Mapping[arrayName] = {
            "firstGhostIndex" : 1,
            "lastGhostIndex" : lastGhostIndex,
            "map" : arrayMapping
        }
    # print(Mapping)

        
    

def FinalizeParticleStruct(*args, **kwargs):
    global FinalStructs, props, arrayLens
    if ndim < 2 or ndim > 3:
        raise RuntimeError("Must set ndim to 2 or 3")

    posProps = [("x", True), ("y", True)]
    props.pop("x")
    props.pop("y")
    if ndim == 3:
        posProps.append(("z", True))
        props.pop("z")
    # TODO: better check to make sure we actually defined the same number and an accidentally defined z doesn't slip in

    real_sol = [(k,v[2]) for k,v in props.items() if v[0] == "real" and v[1] == True and v[4] == False]
    int_sol = [(k,v[2]) for k,v in props.items() if v[0] == "int" and v[1] == True and v[4] == False]

    real_prop = [(k,v[2]) for k,v in props.items() if v[0] == "real" and v[1] == False and v[4] == False]
    int_prop = [(k,v[2]) for k,v in props.items() if v[0] == "int" and v[1] == False and v[4] == False]

    proj = [(k,v[2]) for k,v in props.items() if v[4] == True]
    
    real_sol_wpos = posProps + real_sol

    if len(int_sol) != 0:
        raise RuntimeError("Tried to create an integer solution property")
    # print(real_sol)
    # print(real_sol_wpos)
    # print(int_sol)
    # print(real_prop)
    # print(int_prop)

    def MakeArrayDeclDict(a):
        nonlocal posProps, real_sol, int_sol, real_prop, int_prop, real_sol_wpos, proj
        d = {"name" : a[0]}
    
        if a[1] == "pos":
            if len(posProps) != ndim:
                raise RuntimeError("NDIM doesn't match number xyz defs")
            d["type"] = "real"
            d["n"] = len(posProps)
        elif a[1] == "real_sol":
            d["type"] = "real"
            d["n"] = len(real_sol)
        elif a[1] == "int_sol":
            d["type"] = "int"
            d["n"] = len(int_sol)
        elif a[1] == "real_prop":
            d["type"] = "real"
            d["n"] = len(real_prop)
        elif a[1] == "int_prop":
            d["type"] = "int"
            d["n"] = len(int_prop)
        elif a[1] == "real_sol_wpos":
            d["type"] = "real"
            d["n"] = len(real_sol_wpos)
        elif a[1] == "proj":
            d["type"] = "real"
            d["n"] = len(proj)
        else:
            raise RuntimeError
        
        return d


    if len(real_sol) + len(int_sol) + len(real_prop) + len(int_prop) + len(proj) != len(props.keys()):
        raise RuntimeError("Failed to generate particle struct, particle lengths dont match")

    FinalStructs = {}
    for struct in layout:
        FinalStructs[struct["StructName"]] = [MakeArrayDeclDict(a) for a in struct["arrays"]]


    GenerateMapping(pos=posProps, real_sol=real_sol, int_sol=int_sol, real_prop=real_prop, int_prop=int_prop, real_sol_wpos=real_sol_wpos, proj=proj)
    


def PRINTPROPS():
    for name in props.keys():
        name + ' ' + str(props[name][0]) + " " + str(props[name][1]) + ' ' + str(props[name][2])
        
def ComponentDeclStr(component):
    typestr = lambda t: "real*8" if t == "real" else "integer*4"
    return typestr(component["type"]) + " :: " + component["name"] + "(" + str(component["n"]) + ")"
def GenStructs():
    # print(FinalStructs)
    if FinalStructs is None:
        raise RuntimeError("Must call Finalize before the particle structure is usable.")

    finalCode = ""
    
    for structName, components in FinalStructs.items():
        # print("Generating " + str(structName))
        s = "type :: " + structName + "\n"
        for component in components:
            if component["n"] == 0:
                continue
            s += "    " + ComponentDeclStr(component) + "\n"
        s += "end type " + structName + "\n\n"
        finalCode += s

    return finalCode

def DetermineReferences(**kwargs):
    specifiedStruct = None
    specifiedSubstruct = None

    if "substructures" in kwargs.keys():
        # if type(kwargs["substructures"]) is not list and type(kwargs["substructures"]) is not tuple:
        #     kwargs["substructures"] = [kwargs["substructures"],]
        #     print(type(kwargs["substructures"]) is not list)
        #     raise RuntimeError
        # print(type(kwargs["substructures"]))
        if len(kwargs["substructures"]) == 2:
            specifiedStruct = kwargs["substructures"][0]
            specifiedSubstruct = kwargs["substructures"][1]
        elif len(kwargs["substructures"]) == 1:
            if kwargs["substructures"][0] in [struct["StructName"] for struct in layout]:
                specifiedStruct = kwargs["substructures"][0]
            else:
                specifiedSubstruct = kwargs["substructures"][0]
    else:
        if "specifiedStruct" in kwargs.keys():
            specifiedStruct = kwargs["specifiedStruct"]
        if "specifiedSubstruct" in kwargs.keys():
            specifiedSubstruct = kwargs["specifiedSubstruct"]

    propName = kwargs["propName"]
    
    containingArrays = [k for k, v in Mapping.items() if propName in v["map"].keys()]

    possibleRefs = [(struct, array[0], array[1]) for struct in layout for array in struct["arrays"] if array[1] in containingArrays]
    if specifiedStruct is not None:
        possibleRefs = [r for r in possibleRefs if r[0]["StructName"] == specifiedStruct]
    
    if specifiedSubstruct is not None:
        possibleRefs = [r for r in possibleRefs if specifiedSubstruct in [a[2] for a in r[0]["arrays"] if a[0] == r[1]][0]]

    if len(possibleRefs) == 0:
        # print(Mapping)
        print([(k,v["map"].keys())  for k, v in Mapping.items()])
        print("SpecifiedStruct: " + str(specifiedStruct) + "\nSpecifiedSubStruct: " + str(specifiedSubstruct) + "\nPropName: " + propName)

        raise RuntimeError("Cannot reference particle property " + str(kwargs))
    elif len(possibleRefs) > 1:
        print(possibleRefs)
        print("SpecifiedStruct: " + str(specifiedStruct) + "\nSpecifiedSubStruct: " + str(specifiedSubstruct) + "\nPropName: " + propName)
        raise RuntimeError("Ambiguous particle property reference " + str(kwargs))
    
    return possibleRefs

def GenPropAccess(*args):
    if len(args) < 2:
        raise RuntimeError("Not enough arguments to PropAccess macro")
    particleIndex = args[0]
    propName = args[-1].lower()
    


    

    chosenRef = DetermineReferences(propName=propName, substructures=args[1:-1])[0]

    refStr = chosenRef[0]["ppiclfArrayName"] + "(" + str(particleIndex) + ")%" + chosenRef[1] + "(" + str(Mapping[chosenRef[2]]["map"][propName]) + ")"
    return refStr


def Loop_All_RealArrays():    
    return [(struct["ppiclfArrayName"], array[0], arrayLens[array[1]]) for struct in layout for array in struct["arrays"] if arrayLens[array[1]] != 0 and array[1] in ["pos", "real_sol", "real_sol_wpos", "real_prop", "proj"]]

def Loop_All_IntArrays():    
    return [(struct["ppiclfArrayName"], array[0], arrayLens[array[1]]) for struct in layout for array in struct["arrays"] if arrayLens[array[1]] != 0 and array[1] in ["int_sol", "int_prop"]]

def Loop_All_SlnArrays():
    return [(struct["ppiclfArrayName"], array[0], arrayLens[array[1]]) for struct in layout for array in struct["arrays"] if arrayLens[array[1]] != 0 and array[0] in ["y_real", "pos"]]

def Loop_All_StructArrays():
    return [(struct["StructName"], struct["ppiclfArrayName"]) for struct in layout]

def Loop_All_SLNProps():
    return [p for p in Mapping["real_sol_wpos"]["map"].keys()]


def GenerateArraySlice(*args):
    particleIndex = args[0]
    startProp = args[1].lower()
    endProp = args[2].lower()
    if startProp[0] in "[(":
        startProp = [s.strip() for s in startProp[1:-1].split(",")]
    else:
        startProp = (startProp,)
    if endProp[0] in "[(":
        endProp = [s.strip() for s in endProp[1:-1].split(",")]
    else:
        endProp = (endProp,)

    startPropRef = DetermineReferences(propName=startProp[-1], substructures=startProp[0:-1])[0]
    endPropRef = DetermineReferences(propName=endProp[-1], substructures=endProp[0:-1])[0]

    expectedElements = int(args[3]) if len(args) == 4 else None
    
    # print(startPropRef)
    # print(startProp)
    # print(endPropRef)
    # print(endProp)
    # print(expectedElements)

    if startPropRef[0] != endPropRef[0] or startPropRef[1] != endPropRef[1]:
        raise RuntimeError("Trying to slice between two elements not in the same array")

    startIndex = Mapping[startPropRef[2]]["map"][startProp[-1]]
    endIndex = Mapping[endPropRef[2]]["map"][endProp[-1]]

    if expectedElements is not None and (endIndex - startIndex + 1) != expectedElements:
        raise RuntimeError("Didn't generate a slice with the expected number of elements")

    # return str(startIndex) + ":" + str(endIndex)
    return startPropRef[0]["ppiclfArrayName"] + "(" + str(particleIndex) + ")%" + startPropRef[1] + "(" + str(startIndex) + ":" + str(endIndex) + ")"


def CountRealGhostProps():
    total = 0
    for struct in layout:
        for a in struct["arrays"]:
            if a[1] not in ["pos", "real_sol", "real_sol_wpos", "real_prop", "proj"]:
                continue
            total += Mapping[a[1]]["lastGhostIndex"] - Mapping[a[1]]["firstGhostIndex"] + 1
            # print(a[1] + " " + str(total))
    
    return str(total)


def CountIntGhostProps():
    total = 0
    for struct in layout:
        for a in struct["arrays"]:
            if a[1] not in ["int_sol", "int_prop"]:
                continue
            if a[1] not in Mapping:
                continue
            total += Mapping[a[1]]["lastGhostIndex"] - Mapping[a[1]]["firstGhostIndex"] + 1
            # print(a[1] + " " + str(total))
    
    return str(total)


def Loop_All_Ghost_Real():
    return [(struct["ppiclfArrayName"], array[0], Mapping[array[1]]["firstGhostIndex"], Mapping[array[1]]["lastGhostIndex"], arrayLens[array[1]]) for struct in layout for array in struct["arrays"] if arrayLens[array[1]] != 0 and Mapping[array[1]]["lastGhostIndex"] != 0 and array[1] in ["pos", "real_sol", "real_sol_wpos", "real_prop", "proj"]]


def Loop_All_Ghost_Int():
    return [(struct["ppiclfArrayName"], array[0], Mapping[array[1]]["firstGhostIndex"], Mapping[array[1]]["lastGhostIndex"], arrayLens[array[1]]) for struct in layout for array in struct["arrays"] if arrayLens[array[1]] != 0 and Mapping[array[1]]["lastGhostIndex"] != 0 and array[1] in ["int_sol", "int_prop"]]
