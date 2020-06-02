import pandas as pd
import os
import math
from pyswip import Prolog


def pkt1(dbConn, filename):
    df = pd.read_sql_query('select * from results_PKT1_PDT_uniprot', dbConn)

    with open(filename, "w") as file:
        for index, row in df.iterrows():
            var1 = "pkt_1('{}', '{}', '{}', '{}', {}).".format(row["Perturbagen"], row["Kinase"],
                                                               row["Target"], row["Cell_Line"], row["Score"])
            file.write(var1 + "\n")

def pkt_knowntarget(dbConn, filename):
    df = pd.read_sql_query('select * from KS_relationship where cell_line = "MCF-7" and source="PDT" and confidence > 0', dbConn)

    with open(filename, "w") as file:
        for index, row in df.iterrows():

            var1 = "knowntarget('{}', '{}').".format(row["substrate"], row["kinase"])
            file.write(var1 + "\n")

def klocation(dbConn, filename): # Produces kinase locations
    df = pd.read_sql_query('select Protein, Secretory, Nuclear, Cytosol, Mitochondria from '
                           'Kinase_Loc', dbConn)

    with open(filename, "w") as file:
        for index, row in df.iterrows():
            row1Secretory = float(row["Secretory"].strip('%')) / 100
            row1Nuclear = float(row["Nuclear"].strip('%')) / 100
            row1Cytosol = float(row["Cytosol"].strip('%')) / 100
            row1Mitochondria = float(row["Mitochondria"].strip('%')) / 100

            var1 = "klocation('{}', {}, {}, {}, {}).".format(row["Protein"], row1Secretory,
                                                             row1Nuclear, row1Cytosol, row1Mitochondria)
            file.write(var1 + "\n")


def plocation(dbConn, filename): # Produces substrate locations
    df = pd.read_sql_query('select Protein, Secretory, Nuclear, Cytosol, Mitochondria from '
                           'Protein_Loc', dbConn)

    with open(filename, "w") as file:
        for index, row in df.iterrows():
            var1 = "plocation('{}', {}, {}, {}, {}).".format(row["Protein"], row["Secretory"],
                                                             row["Nuclear"], row["Cytosol"], row["Mitochondria"])
            file.write(var1 + "\n")


def observation(dbConn, filename): # Produces observation probs, p^fc etc.
    df = pd.read_sql_query('select * from Observation where cell_line = "MCF-7"', dbConn)

    # df = df[df["fold_change"] < -1 or df["fold_change"] > 1]

    df = df[df["p_value"] != -888]
    df = df[df["fold_change"] != -888]
    df = df[df["cv"] != -888]

    # max_fold_change = -17.8323360893077

    with open("temp.pl", "w") as file:
        for index, row in df.iterrows():

            # Find absolute log2 fold-change
            absFC = round(abs(row["fold_change"]), 5)

            # Apply formula from (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3957066/)
            # Convert to probability with a p value threshold
            threshold = 0.05
            prob = round(1 - (math.pow(row["p_value"] / threshold, absFC)), 5)
            # Same as putting a 0.05 p value threshold..
            if prob < 0:
                continue

            upOrDown = ''
            if(row["fold_change"] > 0):
                upOrDown = 'up'
            else:
                upOrDown = 'down'

            var1 = "perturbs('{}', '{}', {}, {}).".format(row["perturbagen"], row["substrate"],
                                                          upOrDown, prob)
            file.write(var1 + "\n")



    removeDuplicate("temp.pl", filename)
    os.remove("temp.pl")

def removeDuplicate(input, output):
    # Remove duplicate lines
    lines_seen = set()  # holds lines already seen
    outfile = open(output, "w")
    for line in open(input, "r"):
        if line not in lines_seen:  # not a duplicate
            outfile.write(line)
            lines_seen.add(line)
    outfile.close()

def uniqueness(dbConn, filename): # Produces uniqueness probs of kinase-substrate tuples
    df = pd.read_sql_query(
        'select KS_relationship.substrate, KS_relationship.kinase, Kinase_Loc.Secretory, '
        'Kinase_Loc.Nuclear, Kinase_Loc.Cytosol, Kinase_Loc.Mitochondria '
        'from KS_relationship join Kinase_Loc on Kinase_Loc.Protein = KS_relationship.kinase '
        'where KS_relationship.cell_line = "MCF-7" '
        'order by KS_relationship.substrate', dbConn)

    df.reset_index(inplace=True, drop=True)

    with open(filename, "w") as file:
        for index, row in df.iterrows():

            multiplierList = []
            for index2, row2 in df.iterrows():

                if row2["substrate"] != row["substrate"] or row2["kinase"] == row["kinase"]:
                    continue

                row1Secretory = float(row["Secretory"].strip('%')) / 100
                row1Nuclear = float(row["Nuclear"].strip('%')) / 100
                row1Cytosol = float(row["Cytosol"].strip('%')) / 100
                row1Mitochondria = float(row["Mitochondria"].strip('%')) / 100

                row2Secretory = float(row2["Secretory"].strip('%')) / 100
                row2Nuclear = float(row2["Nuclear"].strip('%')) / 100
                row2Cytosol = float(row2["Cytosol"].strip('%')) / 100
                row2Mitochondria = float(row2["Mitochondria"].strip('%')) / 100

                currMultiplier = (row1Secretory * row2Secretory) + (row1Nuclear * row2Nuclear) + (
                            row1Cytosol * row2Cytosol) + \
                                 (row1Mitochondria * row2Mitochondria)

                multiplierList.append(1 - currMultiplier)

            result = 1
            for mult in multiplierList:
                result *= mult

            var1 = "uniqueness('{}','{}', {}).".format(row["substrate"], row["kinase"], result)
            file.write(var1 + "\n")


def majority(dbConn, filename, recursive = False, consultList = None, iteration = None): # Produces majority probabilities for perturbagen-kinase tuples

    if (recursive):
        getResultsProlog(consultList, iteration, True, dbConn)
        PK = "Return_iteration"
    else:
        PK = "PK_relationship"

    query = 'select {}.kinase, {}.perturbagen, KS_relationship.substrate, observation.fold_change,'\
            'observation.p_value, observation.cv '\
            'from {} join KS_relationship on {}.kinase = KS_relationship.kinase '\
            'join observation on '\
            'observation.perturbagen = {}.perturbagen '\
            'and observation.substrate = KS_relationship.substrate '\
            'and observation.cell_line = KS_relationship.cell_line '\
            'where KS_relationship.cell_line = "MCF-7" '\
            'order by {}.kinase, {}.perturbagen'.format(PK, PK, PK, PK, PK, PK, PK)

    df = pd.read_sql_query(query, dbConn)

    df = df[df["p_value"] != -888]
    df = df[df["fold_change"] != -888]
    df = df[df["cv"] != -888]

    # df = df[df["p_value"] <= 0.05]
    # df = df[~df["fold_change"].between(-1,1)]

    df.reset_index(inplace=True, drop=True)

    with open(filename, "w") as file:
        prevKinase = df["kinase"][0]
        prevPerturbagen = df["perturbagen"][0]
        up = 0
        down = 0

        for index, row in df.iterrows():
            newKinase = row["kinase"]
            newPerturbagen = row["perturbagen"]

            if newKinase == prevKinase and newPerturbagen == prevPerturbagen:
                if row["fold_change"] > 0:
                    up += row["p_value"] ** abs(row["fold_change"])
                elif row["fold_change"] < 0:
                    down += row["p_value"] ** abs(row["fold_change"])
            else:
                if down != 0:
                    majorityRatio = round(up / (down + up), 6)
                else:
                    majorityRatio = 0

                var1 = "majority('{}', '{}', {}).".format(row["perturbagen"], row["kinase"], majorityRatio)
                file.write(var1 + "\n")

                up = 0
                down = 0

                if row["fold_change"] > 0:
                    up += row["p_value"] ** abs(row["fold_change"])
                elif row["fold_change"] < 0:
                    down += row["p_value"] ** abs(row["fold_change"])

                prevKinase = newKinase
                prevPerturbagen = newPerturbagen


def getResultsProlog(consultList, iteration, insertDB = False, dbConn = None):
    prolog = Prolog()

    for file in consultList:
        prolog.consult(file)

    query = ''
    if iteration == 0:
        query = "doesXinhibitAinCzero(Perturbagen, Kinase, Cell_Line, Prob)"
    elif iteration == 1:
        query = "doesXinhibitAinC(Perturbagen, Kinase, Target, Cell_Line, Prob)"
    elif iteration == 2:
        query = "doesXinhibitAinC2(Perturbagen, Kinase, Target, Cell_Line, Prob, PerturbagenProb, ColocProb, UniquenessProb)"
    elif iteration == 4:
        query = "doesXinhibitAinC4(Perturbagen, Kinase, Cell_Line, Prob)"
    elif iteration == 5:
        query = "doesXinhibitAinC5(Perturbagen, Kinase, Cell_Line, Prob)"

    solutions = prolog.query(query)

    solutionlist = []
    for solution in solutions:
        solutionlist.append(solution)

    if(insertDB):
        dbConn.cursor().execute('drop table Return_iteration')

        sql_query_create_table = 'create table Return_iteration (kinase TEXT, perturbagen TEXT, probability TEXT);'
        dbConn.cursor().execute(sql_query_create_table)

        for solution in solutionlist:
            kinase = solution.get('Kinase')
            perturbagen = solution.get('Perturbagen')
            prob = str(solution.get('Prob'))

            sql_query_insert = 'insert into Return_iteration (kinase, perturbagen, probability) values(\"' + kinase + '\",\"' + perturbagen + '\",\"' + prob + '\")'

            dbConn.cursor().execute(sql_query_insert)
            dbConn.commit()

    return solutionlist

def queryProlog(query, consultList):
    prolog = Prolog()

    for file in consultList:
        prolog.consult(file)

    solutions = prolog.query(query)

    solutionlist = []
    for solution in solutions:
        solutionlist.append(solution)


    return solutionlist

def generate_kinase_list(dbConn):

    df = pd.read_sql_query('select distinct kinase from KS_relationship', dbConn)

    with open("kinase.pl", "w") as file:
        for index, row in df.iterrows():
            var1 = "kinase('{}').".format(row["kinase"])
            file.write(var1 + "\n")

def knowninhibitor(dbConn):
    df = pd.read_sql_query('select * from PK_relationship', dbConn)

    with open("knowninhibitor.pl", "w") as file:
        for index, row in df.iterrows():
            var = "knowninhibitor('{}', '{}', '{}','{}')."\
                .format(row["perturbagen"], row["kinase"], row["source"], row["score"])

            file.write(var + "\n")




# # Generates colocation probabilities of each Kinase-Substrate combinations
# # Defunct
# def colocation(dbConn, filename):
#     df = pd.read_sql_query('Select Kinase_Loc.Protein as Kinase, Protein_Loc.Protein as Protein, '
#                            'Kinase_Loc.Secretory as KSec, Kinase_Loc.Nuclear as KNuc, Kinase_Loc.Cytosol as KCyt, '
#                            'Kinase_Loc.Mitochondria as KMit, Protein_Loc.Secretory as PSec, Protein_Loc.Nuclear as PNuc, '
#                            'Protein_Loc.Cytosol as PCyt, Protein_Loc.Mitochondria as PMit '
#                            'from Protein_Loc join Kinase_Loc', dbConn)
#
#     df.reset_index(inplace=True, drop=True)
#
#     with open(filename, "w") as file:
#         for index, row in df.iterrows():
#             ksec, knuc, kcyt, kmit = float(row["KSec"].strip('%')) / 100, float(row["KNuc"].strip('%')) / 100, \
#                                      float(row["KCyt"].strip('%')) / 100, float(row["KMit"].strip('%')) / 100
#
#             psec, pnuc, pcyt, pmit = row["PSec"], row["PNuc"], row["PCyt"], row["PMit"]
#
#             prob = (ksec * psec) + (knuc * pnuc) + (kcyt * pcyt) + (kmit * pmit)
#
#             var1 = "colocation('{}', '{}', {}).".format(row["Protein"], row["Kinase"], prob)
#             file.write(var1 + "\n")

def removeDuplicate(input, output):
    # Remove duplicate lines
    lines_seen = set()  # holds lines already seen
    outfile = open(output, "w")
    for line in open(input, "r"):
        if line not in lines_seen:  # not a duplicate
            outfile.write(line)
            lines_seen.add(line)
    outfile.close()