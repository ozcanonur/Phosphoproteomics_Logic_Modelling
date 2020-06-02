from pyswip import Prolog
import pandas as pd
prolog = Prolog()

# Consult
prolog.consult('Consult/RuleConsultTest.pl')
prolog.consult('Consult/Protein_expression(expressed_in).pl')
prolog.consult('Consult/KS_relationship(known_substrate).pl')
prolog.consult('Consult/Observations_NEW_PROB.pl')
prolog.consult('Consult/plocation.pl')
prolog.consult('Consult/klocation.pl')
prolog.consult('Consult/PonProtein(ponProtein).pl')
prolog.consult('Consult/Uniqueness_NEW_PROB.pl')

prolog.consult('Consult/Majority_NEW_PROB.pl')
# prolog.consult('Consult/Majority_NEW_PROB_3rditerationReturned.pl')

def getSolutions(iteration):
    query = ''
    if iteration == 0:
        query = "doesXinhibitAinCzero(Perturbagen, Kinase, Cell_Line, Prob)"
    elif iteration == 1:
        query = "doesXinhibitAinC(Perturbagen, Kinase, Cell_Line, Prob)"
    elif iteration == 2:
        query = "doesXinhibitAinC2(Perturbagen, Kinase, Cell_Line, Prob)"
    elif iteration == 4:
        query = "doesXinhibitAinC4(Perturbagen, Kinase, Cell_Line, Prob)"
    elif iteration == 5:
        query = "doesXinhibitAinC5(Perturbagen, Kinase, Cell_Line, Prob)"

    solutions = prolog.query(query)

    solutionlist = []
    for solution in solutions:
        solutionlist.append(solution)

    return solutionlist

solutionlist = getSolutions(2)

df = pd.DataFrame(solutionlist)
print(solutionlist)


import sqlite3

conn = sqlite3.connect("chemphopro.db")

def createTableFromIterationResults(conn, solutionlist):
    conn.cursor().execute('drop table Return_iteration')

    sql_query_create_table = 'create table Return_iteration (kinase TEXT, perturbagen TEXT, probability TEXT);'
    conn.cursor().execute(sql_query_create_table)

    for solution in solutionlist:
        kinase = solution.get('Kinase')
        perturbagen = solution.get('Perturbagen')
        prob = str(solution.get('Prob'))

        sql_query_insert = 'insert into Return_iteration (kinase, perturbagen, probability) values(\"' + kinase + '\",\"' + perturbagen + '\",\"' + prob + '\")'

        conn.cursor().execute(sql_query_insert)
        conn.commit()

createTableFromIterationResults(conn, solutionlist)