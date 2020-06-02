import sqlite3
import pandas as pd
import numpy as np
import sys
from datetime import datetime
import itertools
import pickle

def build_obs_dict(perturbagen, p_threshold, cell_line, dict={}):
    conn = sqlite3.connect("omnipath.db")

    queryString = 'select * from Observation ' \
                  'where perturbagen = "{}" and p_value < {} and p_value <> -888.0 and fold_change <> -888.0 ' \
                  'and cell_line = "{}" and substrate not like "%(X%" and substrate not like "%(M%"'.format(perturbagen, p_threshold, cell_line)

    df = pd.read_sql_query(queryString, conn)

    for index, row in df.iterrows():

        entry_info = [row['fold_change'], row['p_value'], row['cv']]
        substrate = row['substrate']
        key = substrate.split('(')[0] + '_' + substrate[substrate.find('(') + 1 : substrate.find(')')]
        dict[key] = entry_info

    return dict


def build_rel_dict(obs_dict, dict={}, include_unknown_and_conflicted=False):

    def remove_unnecessary_single_kinases(dict):
        for key in dict:
            curr_list = dict[key]
            i = 0
            while i < len(curr_list):
                ele = curr_list[i]
                if ele[0][1] != '':
                    i += 1
                    continue
                curr_kinase = ele[0][0]
                k = 0
                while k < len(curr_list):
                    ele2 = curr_list[k]
                    if ele2[0][1] != '' and curr_kinase == ele2[0][0]:
                        curr_list.remove(ele)
                        dict[key] = curr_list
                        i -= 1
                        break
                    k += 1
                i += 1
        return dict

    conn = sqlite3.connect("omnipath.db")

    unknown_and_conflicted = ''
    if not include_unknown_and_conflicted:
        unknown_and_conflicted = 'and PsT_effect <> "unknown" and PsT_effect <> "conflicting" '

    queryString = 'select second_table.KPa as left_prot, second_table.TProtein as mid_prot_2, second_table.PsT as mid_ps_2, ' \
                  'first_table.KPa as mid_prot, first_table.TProtein as right_prot, first_table.PsT as right_ps, ' \
                  'first_table.KPa_class as right_effect, known_sign.PsT_effect ' \
                  'from (select * from known_target where KPa <> TProtein) as first_table left join ' \
                  '(select * from known_target where KPa <> TProtein) as second_table on first_table.KPa = second_table.TProtein ' \
                  'join known_sign on mid_prot = known_sign.KPa and right_ps = known_sign.PsT ' \
                  'where (left_prot is null or left_prot <> right_prot) ' + unknown_and_conflicted + \
                  'group by mid_ps_2, right_ps order by right_ps, mid_ps_2 desc'

    df = pd.read_sql_query(queryString, conn)

    for index, row in df.iterrows():
        right_substrate = row['right_ps']
        right_prot = right_substrate.split('(')[0]
        right_ps = right_substrate[right_substrate.find("(") + 1: right_substrate.find(")")]
        right_substrate = right_prot + '_' + right_ps

        effect = row['right_effect']
        reg = row['PsT_effect']

        mid_substrate = row['mid_ps_2']
        mid_prot = row['mid_prot']

        if right_substrate not in obs_dict:
            continue
        elif (mid_substrate is not None) and (mid_prot + '_' + mid_substrate[mid_substrate.find('(') + 1: mid_substrate.find(')')]) in obs_dict:
            mid_ps = mid_substrate[mid_substrate.find('(') + 1: mid_substrate.find(')')]
            new_entry = [[mid_prot, mid_ps], effect, reg]
        elif mid_substrate is None:
            new_entry = [[mid_prot, ''], effect, reg, 'No_BG']
        else:
            new_entry = [[mid_prot, ''], effect, reg, 'No_EXP']

        if right_substrate in dict:
            if new_entry not in dict[right_substrate]:
                curr_value = dict[right_substrate][:]
                curr_value.append(new_entry)

                dict[right_substrate] = curr_value
        else:
            dict[right_substrate] = [new_entry]

    return remove_unnecessary_single_kinases(dict)


def check_if_exists(node, path, check_ps_in_path=False):
    if check_ps_in_path:
        prot = node[0][0]
        ps = node[0][1]

        for element in path:
            if not isinstance(element, list) and element == prot + '_' + ps:
                return True
            if element[0][0] == prot and element[0][1] == ps:
                return True
        return False
    else:
        prot = node[0][0]

        for element in path:
            if not isinstance(element, list) and element.split('_')[0] == prot:
                return True
            else:
                if element[0][0] == prot:
                    return True
        return False


def makes_sense(start, node, obs_dict):
    node_obs = node[0][0] + '_' + node[0][1]
    if (obs_dict[start][0] < 0 and node[1] == 'K') or (obs_dict[start][0] > 0 and node[1] == 'Pa'):
        return (obs_dict[node_obs][0] > 0 and node[2] == 'p_dec') or (obs_dict[node_obs][0] < 0 and node[2] == 'p_inc')
    elif (obs_dict[start][0] < 0 and node[1] == 'Pa') or (obs_dict[start][0] > 0 and node[1] == 'K'):
        return (obs_dict[node_obs][0] > 0 and node[2] == 'p_inc') or (obs_dict[node_obs][0] < 0 and node[2] == 'p_dec')
    else:
        return False


def add_details_to_node(node, obs_dict, include_obs=False):
    if include_obs:
        info = [[node[0][0], node[0][1]]]
        obs_info = obs_dict[node[0][0] + '_' + node[0][1]][:]
    else:
        info = [[node[0][0], '']]
        obs_info = []

    k_or_pa = node[1]
    reg = node[2]
    info.extend([k_or_pa, reg, obs_info])

    return info


def bottomup_path(rel_dict, obs_dict, start, path=[], check_ps_in_path=False, check_makes_sense=False):

    if isinstance(start, list) and start[0][1] != '':
        info = add_details_to_node(start, obs_dict, include_obs=True)
        start = start[0][0] + '_' + start[0][1]
        path = path + [info]
    elif isinstance(start, list):
        info = add_details_to_node(start, obs_dict)
        info.extend([start[3]])
        start = start[0][0]
        path = path + [info]
    else:
        path = path + [start]

    if start not in rel_dict:
        if len(start.split('_')) == 2:
            last_node = path[len(path) - 1]
            last_node.extend(['Unk/Conf'])
        return [path]

    paths = []
    newpaths = []

    for node in rel_dict[start]:
        if check_makes_sense:
            if not check_if_exists(node, path, check_ps_in_path):
                if node[0][1] == '':
                    newpaths = bottomup_path(rel_dict, obs_dict, node, path, check_ps_in_path, check_makes_sense)
                elif makes_sense(start, node, obs_dict):
                    newpaths = bottomup_path(rel_dict, obs_dict, node, path, check_ps_in_path, check_makes_sense)
            else:
                last_node = path[len(path) - 1]
                last_node.extend(['FB_Loop'])
                newpaths = [path]
            for newpath in newpaths:
                paths.append(newpath)
        else:
            if not check_if_exists(node, path, check_ps_in_path):
                newpaths = bottomup_path(rel_dict, obs_dict, node, path, check_ps_in_path, check_makes_sense)
            else:
                last_node = path[len(path) - 1]
                last_node.extend(['FB_Loop'])
                newpaths = [path]
            for newpath in newpaths:
                paths.append(newpath)

    return paths


# get_bottomup_path('AKT1_S473')
def get_bottomup_path(target):

    # sys.setrecursionlimit(2000)
    print('Obs dict start: ' + str(datetime.now()))
    obs_dict = build_obs_dict(perturbagen='Torin', p_threshold=0.2, cell_line='MCF-7')
    print('Obs dict end: ' + str(datetime.now()))
    print('Rel dict start: ' + str(datetime.now()))
    rel_dict = build_rel_dict(obs_dict, include_unknown_and_conflicted=False)
    print('Rel dict end: ' + str(datetime.now()))

    paths = bottomup_path(rel_dict, obs_dict, target, [], check_ps_in_path=False, check_makes_sense=False)
    paths.sort()
    paths = list(paths for paths, _ in itertools.groupby(paths))

    for path in paths:
        print(path)

    print('Path count: ' + str(len(paths)))
    print('Longest path: ' + check_longest_path(paths))


def check_longest_path(paths):
    max = 0
    for path in paths:
        if len(path) > max:
            max = len(path)
    return str(max)


def simplify_paths(paths):
    new_paths = []

    for path in paths:
        new_path = []
        for node in path:
            if not isinstance(node, list):
                new_path.append([node.split('_')[0], node.split('_')[1]])
            else:
                new_path.append([node[0][0], node[0][1]])
        new_paths.append(new_path)

    return new_paths


def build_paths(obs_dict, rel_dict, target, check_ps_in_path=False, check_makes_sense=False, simple_paths=False):
    paths = bottomup_path(rel_dict, obs_dict, target, [], check_ps_in_path, check_makes_sense)
    paths.sort()
    paths = list(paths for paths, _ in itertools.groupby(paths))

    if simple_paths:
        paths = simplify_paths(paths)

    for path in paths:
        print(path)

    print('Path count: ' + str(len(paths)))
    print('Longest path: ' + check_longest_path(paths))

    return paths


# paths = build_paths(obs_dict, rel_dict, 'AKT1_S473', False, False, simple_paths=False)
# paths = build_paths(obs_dict, rel_dict, '1', False, False, simple_paths=False)


def stop_reasons_dict(obs_dict, rel_dict, dict={}):

    conn = sqlite3.connect("omnipath.db")

    queryString = 'select distinct substrate from Observation where perturbagen = "Torin" and ' \
                  'substrate not like "%(X%" and substrate not like "%(M%" and substrate not like "%HUMAN%"'

    df = pd.read_sql_query(queryString, conn)

    total_no_exp = 0
    total_no_bg = 0
    total_unk_conf = 0
    total_fb_loop = 0

    for index, row in df.iterrows():
        substrate = row['substrate']
        substrate = substrate.split('(')[0] + '_' + substrate[substrate.find('(') + 1: substrate.find(')')]

        if substrate not in rel_dict:
            continue

        print('Start building path for ' + substrate)
        paths = build_paths(obs_dict, rel_dict, substrate, False, False, simple_paths=False)

        no_exp = 0
        no_bg = 0
        unk_conf = 0
        fb_loop = 0
        for path in paths:
            if len(path) == 1:
                continue

            last_node = path[len(path) - 1]
            stop_reason = last_node[len(last_node) - 1]

            if stop_reason == 'No_EXP':
                no_exp += 1
            elif stop_reason == 'No_BG':
                no_bg += 1
            elif stop_reason == 'Unk/Conf':
                unk_conf += 1
            elif stop_reason == 'FB_Loop':
                fb_loop += 1
            else:
                print('error for ' + str(path))

        total_no_exp += no_exp
        total_no_bg += no_bg
        total_unk_conf += unk_conf
        total_fb_loop += fb_loop

        dict[substrate] = [no_exp, no_bg, unk_conf, fb_loop]

    total = total_no_exp + total_no_bg + total_unk_conf + total_fb_loop

    print('Total no_exp: ' + str(total_no_exp) + ', ' + str((total_no_exp / total) * 100) + '%')
    print('Total no_bg: ' + str(total_no_bg) + ', ' + str((total_no_bg / total) * 100) + '%')
    print('Total unk_conf: ' + str(total_unk_conf) + ', ' + str((total_unk_conf / total) * 100) + '%')
    print('Total fb_loop: ' + str(total_fb_loop) + ', ' + str((total_fb_loop / total) * 100) + '%')

    total_list_dict = {}
    total_list_dict['no_exp'] = total_no_exp
    total_list_dict['no_bg'] = total_no_bg
    total_list_dict['unk_conf'] = total_unk_conf
    total_list_dict['fb_loop'] = total_fb_loop

    return dict, total_list_dict


# dict = stop_reasons_dict(obs_dict, rel_dict)

def get_all_perturbagen_stop_reasons(target, cell_line, p_value, result_dict={}):

    conn = sqlite3.connect("omnipath.db")
    queryString = 'select distinct perturbagen from Observation where cell_line = "{}"'.format(cell_line)

    df = pd.read_sql_query(queryString, conn)

    for index, row in df.iterrows():
        perturbagen = row['perturbagen']

        obs_dict = build_obs_dict(perturbagen=perturbagen, p_threshold=p_value, cell_line=cell_line)
        rel_dict = build_rel_dict(obs_dict, include_unknown_and_conflicted=False)

        total_list_dict = stop_reasons_dict(obs_dict, rel_dict)[1]

        result_dict[perturbagen] = total_list_dict

    return result_dict



# all_perturbagens_stop_dict = get_all_perturbagen_stop_reasons('AKT1_S473', 'MCF-7', 0.1)


def get_depth_rankings(paths):
    rankings = {}

    for path in paths:
        i = 0
        while i < len(path):
            node = path[i]
            curr_protein = node[0]

            if curr_protein in rankings:
                rankings[curr_protein] = rankings[curr_protein] + (len(path) - i - 1)
            else:
                rankings[curr_protein] = len(path) - i - 1

            i += 1

    rankings = {k: v for k, v in sorted(rankings.items(), key=lambda item: item[1], reverse=True)}

    numbers_list = []
    for key in rankings:
        numbers_list.append(rankings[key])

    rankings_ext = {}
    i = -1
    prev_max = 0
    for key in rankings:
        if prev_max == rankings[key]:
            rankings_ext[str(i)] = rankings_ext[str(i)] + [key]
        else:
            i += 1
            prev_max = rankings[key]
            rankings_ext[str(i)] = [key]

    rankings_order = {}
    for key in rankings_ext:
        for node in rankings_ext[key]:
            rankings_order[node] = key

    return rankings, rankings_order