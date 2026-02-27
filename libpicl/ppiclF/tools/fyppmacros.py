import copy
#the props dict will have the property name as the key
#and an array with 3 elements, 0 is the datatype, 1 is bool for being a solution, 2 is bool for ghost prop
props = {}
propArrays = {}
propArrayNames = ["real_sol", "real_prop", "int_prop", "proj", "real_sol_ghost", "real_prop_ghost", "int_prop_ghost", "interp", "VU_props"]
FinalStructs = None
Mapping = {}
ndim = -1
definedVars = {}

# define the structures as they would appear to the user if we were using regular fortran derived types

"""
Format for the layout specification
{
    # user facing structure name, this is how it will be refered to in macros. Can conflict with other structures in code, but probably shouldn't for readability
    "UserStructName": "PPICLF_t_particle",
    # components of the structure from the user perspective. Key is the name that will refer to it, and the value is the property type that will be put in the array.
    "UserComponents":
    {
        "y":        "real_sol",
        "ydot":     "real_sol",
        "ydotc":    "real_sol",
        "y1":       "real_sol",
        "rprop":    "real_prop",
        "iprop":    "int_prop",
        "feedback": "proj"
    }
    # FortranStructs defines a list of the structs that will exist from the perspective of fortran/the compiler,
    # as well as how the user components defined above will map into them.
    "FortranStructs":
    [
        {
            # this is the actual name of the stucture that will be in the code, this  
            "FortStructName": "PPICLF_t_particlepos",
            "components":
            [
                
            ]
        },
        {
            "FortStructName": "PPICLF_t_particle"
        }
    ]
}
"""
userLayout = [
    {
        "UserStructName": "PPICLF_t_particle",
        "UserComponents":
        {
            "y":        "real_sol",
            "ydot":     "real_sol",
            "ydotc":    "real_sol",
            "y1":       "real_sol",
            "rprop":    "real_prop",
            "iprop":    "int_prop",
            "feedback": "proj"
        },
        "FortranStructs":
        {
            "PPICLF_t_particlepos" : {
                "components":
                {
                    "y":        "PROP,pos"
                }
            },
            "PPICLF_t_particle_data" : {
                "components":
                {
                    "y":        "remaining",
                    "ydot":     "ALL",
                    "ydotc":    "ALL",
                    "y1":       "ALL",
                    "rprop":    "ALL",
                    "iprop":    "ALL",
                    "feedback": "ALL"
                }
            }
        }
    },
    {
        "UserStructName": "PPICLF_t_ghostParticle",
        "UserComponents":
        {
            "y":        "real_sol_ghost",
            "rprop":    "real_prop_ghost",
            "iprop":    "int_prop_ghost",
        },
        "FortranStructs":
        {
            "PPICLF_t_particle_ghost" : {
                "components":
                {
                    "y":        "ALL",
                    "rprop":    "ALL",
                    "iprop":    "ALL"
                }
            }
        }
    }
]

def EvalPropArgs(args):
    solution = False
    savedProp = False
    ghost = False
    vector = False
    interp = False
    proj = False
    dtype = None
    for arg in args:
        arg = arg.lower()
        if arg in ["solution"]:
            solution = True
            continue
        if arg in ["ghost"]:
            ghost = True
            continue
        if arg in ["prop", "saved"]:
            savedProp = True
            continue
        if arg in ["interp", "interpolated"]:
            interp = True
            continue
        if arg in ["proj"]:
            proj = True
            dtype = "real"
            continue
        typeinfo = arg.split("_")
        # print(typeinfo)
        if "vec" in typeinfo:
            vector = True
        if "real" in typeinfo:
            dtype = "real"
        if "int" in typeinfo:
            dtype = "int"
        

    
    if dtype is None:
        raise RuntimeError("Didn't specify type of property")
    
    return {
        "sol": solution,
        "savedProp": savedProp,
        "ghost": ghost,
        "vector": vector,
        "interp": interp,
        "proj": proj,
        "dtype": dtype
    }


def convertStringListToList(string):
    if string[0] in "([":
        string = [s.strip() for s in string[1:-1].split(",")]
    else:
        string = [string]
    return string

def setvar(name, value):
    globals()[name] = value


def AddProp(name, *args):
    name = convertStringListToList(name.lower())
    componentNames = None
    if len(name) > 1:
        componentNames = name[1:]
    name = name[0]
    if name in props.keys():
        raise RuntimeError("Property " + name + " redefined")

    newProp = {"name" : name, "componentNames": componentNames,  **EvalPropArgs(args)}
    # print(newProp)
    props.update({name: newProp})


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

def SelectFortStruct(prop, fortStructs):
    for i in range(len(fortStructs)):
        if type(fortStructs[i][1]) is not list and type(fortStructs[i][1]) is not tuple:
            fortStructs[i][1] = [fortStructs[i][1],]
    for struct in fortStructs:
        correctStruct = True
        for modifier in struct[1]:
            
            if modifier.lower() == "all" or modifier.lower() == "remaining":
                return struct[0]
            if modifier[0:10].lower() == "less,prop," and modifier[10:].lower() == prop["name"]:
                correctStruct = False
                break
            if modifier[0:5].lower() =="prop," and modifier[5:].lower() == prop["name"]:
                return struct[0]
    raise RuntimeError("failed to find an array to put " + str(prop["name"]) + " in")

def DumpFinalStructs():
    global FinalStructs
    print("Final Structs:")
    for mainStructName, mainStruct in FinalStructs.items():
        print(mainStructName)
        # print(mainStruct["components"])
        print("\tExternal Component Arrays")
        for arrayName, includedProps in mainStruct["components"].items():
            print("\t\t" + arrayName + " " + str([pName + str("(" + ",".join(mainStruct["FortranStructs"][p[0]]["components"][arrayName][pName]["componentNames"]) + ")" if mainStruct["FortranStructs"][p[0]]["components"][arrayName][pName]["vector"] else "") for pName, p in includedProps.items()]))
        print("\tInternal Fortran Structures: ")
        for fortStructName, fortStruct in mainStruct["FortranStructs"].items():
            print("\t\t" + fortStructName)
            # print(fortStruct)
            for arrayName, includedProps in fortStruct["components"].items():
                # print(includedProps)
                print("\t\t\t" + arrayName + " " + str([p + "(" + str(v["startIndex"]) + ":" + str(v["startIndex"] + v["len"] - 1) + ")" for p, v in includedProps.items() if type(p) is not int]) + "Total Len: " + str(includedProps[0]))


def FinalizeParticleStruct(*args, **kwargs):
    global FinalStructs, props, propArrays, propArrayNames
    if ndim < 2 or ndim > 3:
        raise RuntimeError("Must set ndim to 2 or 3")
    
    for basePropArray in propArrayNames:
        propArrays[basePropArray] = []
    

    # print(props)
    for prop in props.values():
        if prop["vector"]:
            if prop["componentNames"] is None:
                prop["componentNames"] = ["x", "y", "z"][0:ndim]
            prop["len"] = len(prop["componentNames"])
        else:
            prop["len"] = 1


        used = False
        if prop["sol"] and prop["dtype"] == "real":
            propArrays["real_sol"].append(prop)
            if prop["ghost"]:
                propArrays["real_sol_ghost"].append(prop)
            used=True
        if prop["savedProp"] and prop["dtype"] == "real":
            propArrays["real_prop"].append(prop)
            if prop["ghost"]:
                propArrays["real_prop_ghost"].append(prop)
            used=True
        if prop["savedProp"] and prop["dtype"] == "int":
            propArrays["int_prop"].append(prop)
            if prop["ghost"]:
                propArrays["int_prop_ghost"].append(prop)
            used=True
        if prop["interp"]:
            propArrays["interp"].append(prop)
            used=True
        if prop["proj"]:
            propArrays["proj"].append(prop)
            used=True


        if not used:
            raise RuntimeError("Property " + prop["name"] + " not allocated to any array. " + str(prop))
    # for k,v in propArrays.items():
    #     print(k + ": " + str([ p["name"] for p in v]))

    FinalStructs = {}
    for UserStruct in userLayout:
        # print(UserStruct["UserStructName"])
        final = {}
        final["UserName"] = UserStruct["UserStructName"]
        final["components"] = {}
        final["FortranStructs"] = { fsName : {"components" : {a : {} for a in fs["components"].keys()}} for fsName, fs in UserStruct["FortranStructs"].items()}
        # print(final["FortranStructs"])
        for name, array in UserStruct["UserComponents"].items():
            final["components"][name] = {}
            baseArray = propArrays[array]
            fortArrays = []
            for fortStructName, fortStruct in UserStruct["FortranStructs"].items():
                if name in fortStruct["components"].keys():
                    fortArrays.append([fortStructName, fortStruct["components"][name]])

            for prop in baseArray:
                selectedStruct = SelectFortStruct(prop, fortArrays)
                final["FortranStructs"][selectedStruct]["components"][name].update({prop["name"] : copy.deepcopy(prop)})
                final["components"][name][prop["name"]] = (selectedStruct, len(final["FortranStructs"][selectedStruct]["components"][name]) - 1)
        
        for fortStructName, fortStruct in final["FortranStructs"].items():
            # print("\t" + fortStructName)
            for arrayName, arrayProps in fortStruct["components"].items():
                if arrayName == "rprop3":
                    print(arrayProps)
                    print(arrayProps.items())
                i = 1
                # print(arrayProps)
                for p in arrayProps.values():
                    p["startIndex"] = i
                    i += p["len"]
                # store the array's length in index 0, since all of the properties are stored with string valued keys this wont be mistaken for a property name
                arrayProps[0] = i - 1
                # store the array's type in index 1, for same reasons
                if len(arrayProps) > 1:
                    arrayProps[1] = [v for p, v in arrayProps.items() if p != 0][0]["dtype"]
                else:
                    # just specify a dytpe, it doesn't matter what as the array will be empty anyways
                    arrayProps[1] = "real"

                # print(arrayProps[1])
                # print("\t\t\t", end="")
                # for p in arrayProps.values():
                #     print(p["name"] + "(" + str(p["startIndex"]) + ":" + str(p["startIndex"] + p["len"] - 1) + ")", end=" ")
                # print("")

        FinalStructs[UserStruct["UserStructName"]] = final

    # DumpFinalStructs()


def PRINTPROPS():
    for name in props.keys():
        name + ' ' + str(props[name][0]) + " " + str(props[name][1]) + ' ' + str(props[name][2])
        
def ComponentDeclStr(name, component):
    typestr = lambda t: "real*8" if t == "real" else "integer*4"
    return typestr(component[1]) + " :: " + name + "(" + str(component[0]) + ")"

def GenStructs():
    # return
    # print(FinalStructs)
    # DumpFinalStructs()
    if FinalStructs is None:
        raise RuntimeError("Must call Finalize before the particle structure is usable.")

    finalCode = ""
    
    for structName, arrays in [(fortStructName, fortStruct["components"]) for UserStruct in FinalStructs.values() for fortStructName, fortStruct in UserStruct["FortranStructs"].items()]:
        s = "type :: " + structName + "\n"
        for arrayName, array in arrays.items():

            if array[0] == 0:
                continue
            s += "    " + ComponentDeclStr(arrayName, array) + "\n"
        s += "end type " + structName + "\n\n"
        finalCode += s

    return finalCode

def ListComponents(*args):
    index = None
   
    if len(args) == 2:
        structType = args[0]
        name = args[1]
    elif len(args) == 1 and (args[0] in FinalStructs.keys()):
        return list(FinalStructs[args[0]]["FortranStructs"].keys())
    else:
        # print(args)
        reference = DetermineReferences(args[0])
        name = reference[0][0]
        if len(reference[0]) > 1:
            index = reference[0][1]
        structType = definedVars[name]
    components = []
    for fortStructName in FinalStructs[structType]["FortranStructs"].keys():
        components.append(name + "__" + fortStructName + (index if index is not None else ""))
    return components

def DeclarePartVar(structType, name, *args):
    global definedVars
    definedVars[name] = structType
    finalCode = ""
    for fortStructName in FinalStructs[structType]["FortranStructs"].keys():
        finalCode += "type(" + fortStructName + ") :: " + name + "__" + fortStructName + "".join(args) + "\n"
    
    return finalCode


def ImportModVar(structType, name):
    global definedVars
    definedVars[name] = structType
    return ", ".join(ListComponents(structType, name))

def DetermineReferences(refString):
    components = refString.split("%")
    # print(refString)
    # print(components)
    for i, component in enumerate(components):
        if "(" not in component:
            components[i] = [component]
            continue
        idx = component.index("(")
        name = component[0:idx]
        index = component[idx:]
        components[i] = [name, index]
            
    if components[0][0] not in definedVars:
        raise ValueError(components[0][0] + " is not defined.")
    return components
    

def GenPropAccess(givenRef, *args):
    # print("got here " + givenRef)
    UserReference = DetermineReferences(givenRef)
    # print(givenRef)
    # print(UserReference)
    try:
        UserStruct = definedVars[UserReference[0][0]]
    except Exception as e:
        print("Lookup of " + str(UserReference[0][0] + " failed. " + givenRef))
    # hacky way to access a whole array, only works if every element in the user array is stored in the same fortran struct
    skipPropIndex = False
    if len(UserReference) == 2:
        anyPropName = list(FinalStructs[UserStruct]["components"][UserReference[1][0]].keys())[0]
        UserReference.append([anyPropName, ])
        skipPropIndex = True
    
    if len(args) != 0:
        if "skipIndex" in args:
            skipPropIndex = True
        else:
            raise ValueError("invalid arguments to USEPARTICLE: " + str(args))
    UserReference[2][0] = UserReference[2][0].lower()
    fortStructName = FinalStructs[UserStruct]["components"][UserReference[1][0]][UserReference[2][0]][0]

    prop = FinalStructs[UserStruct]["FortranStructs"][fortStructName]["components"][UserReference[1][0]][UserReference[2][0]]

    fortRef = UserReference[0][0] + "__" + fortStructName # mangled name of variable
    fortRef += UserReference[0][1] if len(UserReference[0]) > 1 else ""
    fortRef += "%" + UserReference[1][0] # component array
    # fortRef += "%" + UserReference[2][0]
    if skipPropIndex:
        pass
    elif prop["vector"] and len(UserReference) > 3:
        fortRef += "(" + str(prop["startIndex"] + prop["componentNames"].index(UserReference[3][0].lower())) + ")"
    elif prop["vector"]: # and len(UserReference) == 3: (always true after first if)
        fortRef += "(" + str(prop["startIndex"]) + ":" + str(prop["startIndex"] + prop["len"] -1 )+ ")"
    else: # non vector component
        fortRef += "(" + str(prop["startIndex"]) + ")"
    # print(fortRef)
    return fortRef
    
def Loop_All_Reals(*givenRefs, sameArrays=False):
    def GenRef(UserRefs, userStructs, fortStructs, arrayNames, endPropName, startIndexes):
        endIndexes = [ FinalStructs[userStructs[i]]["FortranStructs"][fortStructs[i]]["components"][arrayNames[i]][endPropName]["startIndex"] + FinalStructs[UserStructs[i]]["FortranStructs"][fortStructs[i]]["components"][arrayNames[i]][endPropName]["len"] -1 for i in range(len(givenRefs))]
        runLen = endIndexes[0] - startIndexes[0] + 1
        # print("swapped runlen(" +  str(runLen) + ") " + str(list(zip(swappedStructs, startIndexes, endIndexes))))
        ref = [runLen]
        for i in range(len(givenRefs)):
            refStr = UserRefs[i][0][0] + "__" + fortStructs[i] + (UserRefs[i][0][1] if len(UserRefs[i][0]) > 1 else "")
            refStr += "%" + arrayNames[i]
            ref += [refStr, startIndexes[i] - 1]
        # print("added ref " + str(ref))
        return ref

    if len(givenRefs) == 1:
        return Loop_All_Reals_SingleRef(givenRefs[0])
    # print(givenRefs)
    UserRefs = [DetermineReferences(given) for given in givenRefs]
    UserStructs = [definedVars[ur[0][0]] for ur in UserRefs ]
    if sameArrays:
        sharedArrayNames = [[arrayName] * len(givenRefs) for arrayName in FinalStructs[UserStructs[0]]["components"].keys() if all([arrayName in list(FinalStructs[us]["components"].keys()) for us in UserStructs])]
    else:
        sharedArrayNames = [[UserRefs[i][1][0] for i in range(len(givenRefs))]]
    print(sharedArrayNames)
    refs = []
    for arrayNames in sharedArrayNames:
        if [fs["components"][arrayNames[0]][1] for fsName, fs in FinalStructs[UserStructs[0]]["FortranStructs"].items() if arrayNames[0] in list(fs["components"].keys())][0] != "real":
            continue
        startStructs = [None] * len(givenRefs)
        startIndexes = None
        prevPropName = None
        jMax = len(FinalStructs[UserStructs[0]]["components"][arrayNames[0]].keys()) - 1
        for j, propName in enumerate(FinalStructs[UserStructs[0]]["components"][arrayNames[0]].keys()):
            sharedProp = all([propName in list(FinalStructs[UserStructs[i]]["components"][arrayNames[i]].keys()) for i in range(len(givenRefs))])
            if startIndexes is None and not sharedProp:
                # print("full skipping " + str(arrayNames) + "%" + propName)
                continue
                # print("\t" + str([propName in list(FinalStructs[UserStructs[i]]["components"][arrayName].keys()) for i in range(len(givenRefs))]))
            
            # end current reference set without starting a new one, as at least one struct doesn't contain the current property
            if not sharedProp:
                # print("ending run. last included prop " + str(prevPropName) + " excluded starts at " + propName)
                refs.append(GenRef(UserRefs, UserStructs, startStructs, arrayNames, prevPropName, startIndexes))
                startIndexes = None
                continue

            fortStructs = [FinalStructs[UserStructs[i]]["components"][arrayNames[i]][propName][0] for i in range(len(givenRefs))]
            if startIndexes is None:
                startIndexes = [ FinalStructs[UserStructs[i]]["FortranStructs"][fortStructs[i]]["components"][arrayNames[i]][propName]["startIndex"] for i in range(len(givenRefs))]
                startStructs = fortStructs
            
            # print(str(j) + " " + propName)
            # print(list(zip(startStructs, startIndexes)))
            
            swappedStructs = [startStructs[i] != fortStructs[i] for i in range(len(givenRefs))]
            # end current reference set and start a new one, because one of them changed internal structure
            if any(swappedStructs):
                refs.append(GenRef(UserRefs, UserStructs, startStructs, arrayNames, prevPropName, startIndexes))
                startIndexes = [ FinalStructs[UserStructs[i]]["FortranStructs"][fortStructs[i]]["components"][arrayNames[i]][propName]["startIndex"] for i in range(len(givenRefs))]


            prevPropName = propName
            startStructs = fortStructs
        
        # add the last ref
        if startIndexes is not None:
            refs.append(GenRef(UserRefs, UserStructs, startStructs, arrayNames, prevPropName, startIndexes))

    # for r in refs:
    #     print(r)

    return refs

def Loop_All_Reals_SingleRef(givenRef):
    UserReference = DetermineReferences(givenRef)
    UserStruct = definedVars[UserReference[0][0]]
    realArrays = []
    # print(UserReference)
    for fortStructName, fortStruct in FinalStructs[UserStruct]["FortranStructs"].items():
        for arrayName, array in fortStruct["components"].items():
            if array[1] != "real": # array's data type
                continue
            if len(UserReference) > 1 and arrayName != UserReference[1][0]:
                continue
            realArrays.append((fortStructName, arrayName, array[0]))
    # print(realArrays)
    retVals = []
    for fortStructName, arrayName, n in realArrays:
        fortRefs = []
        fortRefs.append((UserReference[0][0] + "__" + fortStructName + (UserReference[0][1] if len(UserReference[0]) > 1 else "")))
        fortRefs.append(arrayName)
        retVals.append(("%".join(fortRefs), n))
    return retVals

def Loop_All_Ints(givenRef):
    UserReference = DetermineReferences(givenRef)
    UserStruct = definedVars[UserReference[0][0]]
    realArrays = []
    for fortStructName, fortStruct in FinalStructs[UserStruct]["FortranStructs"].items():
        for arrayName, array in fortStruct["components"].items():
            if array[1] != "int": # array's data type
                continue
            realArrays.append((fortStructName, arrayName, array[0]))
    # print(realArrays)
    retVals = []
    for fortStructName, arrayName, n in realArrays:
        fortRefs = []
        fortRefs.append((UserReference[0][0] + "__" + fortStructName + (UserReference[0][1] if len(UserReference[0]) > 1 else "")))
        fortRefs.append(arrayName)
        retVals.append(("%".join(fortRefs), n))
    return retVals

def Loop_All_Real_Overlaps(givenRef, overlapStruct):
    UserReference = DetermineReferences(givenRef)
    UserStruct = definedVars[UserReference[0][0]]
    refs = []
    for arrayName, array in FinalStructs[UserStruct]["components"].items():
        sharedProps = []
        if arrayName not in FinalStructs[overlapStruct]["components"].keys():
            continue
        if FinalStructs[UserStruct]["FortranStructs"][list(array.values())[0][0]]["components"][arrayName][1] != "real":
            continue
        for propName, propLoc in array.items():
            if propName not in FinalStructs[overlapStruct]["components"][arrayName].keys():
                continue
            sharedProps.append((propName, propLoc))

        skip = None
        sharedPropsFortStructs = {fortStructName : [p for p in sharedProps if p[1][0] == fortStructName] for fortStructName in FinalStructs[UserStruct]["FortranStructs"].keys()}
        # print(sharedPropsFortStructs)
        for fortStructName, sharedProps in sharedPropsFortStructs.items():
            # print(sharedProps)
            i = 0
            startIndex = None
            startProp = None
            currLen = 0
            # print(len(sharedProps))
            while i < len(sharedProps):
                if startIndex is None:
                    startIndex = FinalStructs[UserStruct]["FortranStructs"][fortStructName]["components"][arrayName][sharedProps[i][0]]["startIndex"]
                    startProp = sharedProps[i][0]
                    # print("started " + startProp + " " + str(startIndex))

                currLen += FinalStructs[UserStruct]["FortranStructs"][fortStructName]["components"][arrayName][sharedProps[i][0]]["len"]

                # print("continued " + startProp + " " + str(startIndex) + " " + sharedProps[i][0] + " " + str(currLen))

                if i + 1 == len(sharedProps) or sharedProps[i][1][1] + 1 != sharedProps[i+1][1][1]:
                    # print("ended " + startProp + " " + str(startIndex) + " " + sharedProps[i][0] + " " + str(currLen))
                    refStr = (UserReference[0][0] + "__" + fortStructName + (UserReference[0][1] if len(UserReference[0]) > 1 else ""))
                    refStr += "%" + arrayName + "(" + str(startIndex) + ":" + str(startIndex + currLen - 1) + ")"
                    refs.append((refStr, currLen))
                    # print(refStr)
                    if i + 1 == len(sharedProps):
                        break
                    startIndex = None
                    currLen = 0
                i += 1

    return refs


def CountReals(structType):
    total = 0
    for fortStructName, fortStruct in FinalStructs[structType]["FortranStructs"].items():
        for arrayName, array in fortStruct["components"].items():
            if array[1] == "real":
                total += array[0]
    
    return str(total)

def CountInts(structType):
    total = 0
    for fortStructName, fortStruct in FinalStructs[structType]["FortranStructs"].items():
        for arrayName, array in fortStruct["components"].items():
            if array[1] == "int":
                total += array[0]
    
    return str(total)
