"""
Modified from repo https://github.com/moyiming1/Retrosynthesis-pathway-ranking
  * Increase timeout to 10 s
  * Don't import joblib routines
  * Comment-out if-statement at the bottom
  * Formated with black
"""
from rdkit.Chem import Descriptors
from rdkit import Chem
import signal

# from joblib import Parallel, delayed


class timeout:
    """
    Function for counting time. If a process takes too long to finish, it will be exited by this function.
    """

    def __init__(self, seconds=1, error_message="Timeout"):
        self.seconds = seconds

        self.error_message = error_message

    def handle_timeout(self, signum, frame):
        raise TimeoutError(self.error_message)

    def __enter__(self):
        signal.signal(signal.SIGALRM, self.handle_timeout)
        signal.alarm(self.seconds)

    def __exit__(self, type, value, traceback):
        signal.alarm(0)


def find_largest_mol(chemical_list):
    """
    Find the largest molecule in the supplied chemcial list
    """
    # supply the chemical list
    MW = 0
    large_chemical = ""
    large_ID = 0
    for chemical in chemical_list:
        mol = Chem.MolFromSmiles(chemical["smiles"])
        if Descriptors.ExactMolWt(mol) > MW:
            large_chemical = chemical["smiles"]
            large_ID = chemical["ID"]
            MW = Chem.Descriptors.ExactMolWt(mol)
    out = {"smiles": large_chemical, "ID": large_ID}
    return out


class syn_tree:
    """
    The pathway tree object
    """

    def __init__(self, chemical, reactions, level=0):
        self.is_leaf = False
        self.level = level
        self.smiles = chemical["smiles"]
        self.ID = chemical["ID"]
        self.record_ID = (
            -1
        )  # it is the reaction ID, used to check there is cyclic tree inside
        self.record_data = None
        self.child = []

        if self.level < 25:
            self.expand(reactions)

    def expand(self, reactions):
        results = []
        for reaction in reactions:
            largest_prd = find_largest_mol(reaction["Products"])
            if self.ID == largest_prd["ID"]:
                results.append(reaction)
        if not results:
            self.is_leaf = True
        else:
            # it is a little problematic since we only use the first one
            self.record_data = results[0]["data"]
            self.record_ID = results[0]["ReactionID"]
            self.child = []
            for rct in results[0]["Reactants"]:
                self.child.append(syn_tree(rct, reactions, level=self.level + 1))

    def output_tree(self):
        output = {
            "smiles": self.smiles,
            "ID": self.ID,
            "record_data": self.record_data,
            "record_ID": self.record_ID,
            "child": [c.output_tree() for c in self.child],
        }
        return output

    def find_depth(self):
        return max(self.depth())

    def depth(self, level=0):
        depth = [level]
        for child in self.child:
            depth.extend(child.depth(level + 1))
        return depth

    def check_cyclic(self):
        ID_list = self.get_reactionID()
        if len(ID_list) == len(set(ID_list)):
            return True
        else:
            print(ID_list)
            return False

    def get_reactionID(self):
        ID_list = []
        if self.record_ID != -1:
            ID_list.append(self.record_ID)
        for child in self.child:
            ID_list.extend(child.get_reactionID())
        return ID_list


def build_lib(lib, freq, new_item, role="rct"):
    if new_item in lib:
        index = lib.index(new_item)
        if role == "rct":
            freq[index][0] += 1
        else:
            # there is an error for previous data processing where freq[index][0] += 1 was used
            # but it seems there is no effect on the final trees output
            freq[index][1] += 1

    else:
        lib.append(new_item)
        index = len(lib) - 1
        if role == "rct":
            freq.append([1, 0])
        else:
            freq.append([0, 1])
    return lib, freq, index


def canonicalize(smi):
    mol = Chem.MolFromSmiles(smi)
    for atom in mol.GetAtoms():
        atom.ClearProp("molAtomMapNumber")
    return Chem.MolToSmiles(mol)


def extract_one_patent(record, patentID):
    lib = []
    freq = []
    # print(patentID)
    reactions = []

    for i, line in enumerate(record):
        try:
            rct, rea, prd = line["data"]["smiles"].split(">")
            rcts = canonicalize(rct).split(".")
            prds = canonicalize(prd).split(".")
        except:
            continue

        reactant = []
        product = []
        for item in rcts:
            lib, freq, index = build_lib(lib, freq, item, role="rct")
            reactant.append({"smiles": item, "ID": index})

        for item in prds:
            lib, freq, index = build_lib(lib, freq, item, role="prd")
            product.append({"smiles": item, "ID": index})

        reactions.append(
            {"ReactionID": i, "Reactants": reactant, "Products": product, "data": line}
        )

    # find roots
    roots = []
    for i, row in enumerate(freq):
        if row[0] == 0 and row[1] == 1:
            roots.append(i)

    # build synthesis tree
    trees = []
    # print(len(roots))
    if len(roots) > 100:
        print(
            "Number of roots: {} in patent {}.".format(
                len(roots), patentID #.split(".")[0]
            )
        )
    for root in roots:
        try:
            with timeout(seconds=10):
                root_node = {"smiles": lib[root], "ID": root}
                tree = syn_tree(root_node, reactions)
                if tree.check_cyclic():
                    if tree.find_depth() > 1:
                        trees.append(
                            {"tree": tree.output_tree(), "depth": tree.find_depth()}
                        )
                else:
                    print("Cyclic pathway found. Drop this cyclic pathway.")
        except Exception as e:
            print("error is {}, patent ID is {}".format(e, patentID))
            pass
    return {"patentID": patentID, "trees": trees}


"""
if __name__ == "__main__":
    import pickle

    with open("patent_extraction_testing_data.pkl", "rb") as f:
        data = pickle.load(f)
    extracted_pathways = Parallel(n_jobs=-1, verbose=1)(
        delayed(extract_one_patent)(data[key], key) for key in data.keys()
    )
"""
