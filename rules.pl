%-----------------------------------------------------------------------------------------------------------------------------------%
% First iteration of Perturbagen > Kinase > Target association
:- table pkt_1/5.
pkt_1(Perturbagen, Kinase, Target, Cell_Line, Score) :- 
   expressedin(Kinase, Cell_Line), % KINASE iS expressed in at least one cell line
   perturbs(Perturbagen, Target, down, Perturbs_Score), % A PERTURBAGEN inhibits a TARGET with a score of
   knowntarget(Target, Kinase), % TARGET is associated with (target of) KINASE
   colocalised(Kinase, Target, Colocalisation_Score), % TARGET and KINASE are colocalised with a score of
   Score is Perturbs_Score * Colocalisation_Score.

% Dynamic with asserta
:- dynamic pkt_3/5.
pkt_1_dynamic(Perturbagen, Kinase, Target, Cell_Line, Score) :- 
   expressedin(Kinase, Cell_Line), % KINASE iS expressed in at least one cell line
   perturbs(Perturbagen, Target, down, Perturbs_Score), % A PERTURBAGEN inhibits a TARGET with a score of
   knowntarget(Target, Kinase), % TARGET is associated with (target of) KINASE
   colocalised(Kinase, Target, Colocalisation_Score), % TARGET and KINASE are colocalised with a score of
   Score is Perturbs_Score * Colocalisation_Score, % PERTURBAGEN inhibits KINASE that affects TARGET with a score of
   asserta(pkt_1(Perturbagen, Kinase, Target, Cell_Line, Score)).

% Naive
pkt_1_naive(Perturbagen, Kinase, Target, Cell_Line, Score) :- 
   expressedin(Kinase, Cell_Line), % KINASE iS expressed in at least one cell line
   perturbs(Perturbagen, Target, down, Perturbs_Score), % A PERTURBAGEN inhibits a TARGET with a score of
   knowntarget(Target, Kinase), % TARGET is associated with (target of) KINASE
   colocalised(Kinase, Target, Colocalisation_Score), % TARGET and KINASE are colocalised with a score of
   Score is Perturbs_Score * Colocalisation_Score.

% Probability of the kinase and the target being at the same cellular location
colocalised(Kinase, Target, Score) :-
	ponProtein(Target, Protein),
	plocation(Protein, P_Secretory, P_Nuclear, P_Cytosol, P_Mitochondria), % Target location
	plocation(Kinase, K_Secretory, K_Nuclear, K_Cytosol, K_Mitochondria), % Kinase location
	Score is (P_Secretory * K_Secretory) + (P_Nuclear * K_Nuclear) + (P_Cytosol * K_Cytosol) + (P_Mitochondria * K_Mitochondria). % Score is the probability of the kinase and the target being at the same cellular location
%-----------------------------------------------------------------------------------------------------------------------------------%
% Second iteration of Perturbagen > Kinase > Target association
pkt_2(Perturbagen, Kinase, Target, Effect, Cell_Line, Perturbs_Score, Coloc_Score, KT_Assoc_Score, Score) :- 
   expressedin(Kinase, Cell_Line), % KINASE is expressed in at least one cell line
   perturbs(Perturbagen, Target, Effect, Perturbs_Score), % A PERTURBAGEN affects a TARGET with a score of
   knowntarget(Target, Kinase), % TARGET is associated with (target of) KINASE
   colocalised(Kinase, Target, Coloc_Score), % TARGET and KINASE are colocalised with a score of
   kt_assoc(Kinase, Target, Cell_Line, KT_Assoc_Score, Single_K_Score), % KINASE is uniquely affecting the TARGET with a score of
   Score is (Perturbs_Score * Coloc_Score) * KT_Assoc_Score. % PERTURBAGEN inhibits KINASE that affects TARGET with a score of

% KT_Assoc ¯\_(ツ)_/¯
kt_assoc(Kinase, Target, Cell_Line, Score, Single_K_Score) :-
   k_affects_t_with_avg_score(Kinase, Target, Cell_Line, Single_K_Score), % Average inhibition score of all PERTURBAGENS affecting a KINASE on a TARGET	
   findall(Kinase1, knowntarget(Target, Kinase1), KinaseList), % KinaseList is the list of all valid KINASES that are potentially affecting a TARGET
   all_k_that_affect_t_with_sum_of_avg_score(KinaseList, Kinase, Target, Cell_Line, Sum), % Sum is the sum of all PERTURBAGENS affecting all valid KINASES that are potentially affecting a TARGET
   (Single_K_Score is 0 -> Score is 0 ; Score is Single_K_Score / Sum). % Score is the ratio of average inhibition score of all PERTURBAGENS affecting a KINASE on a TARGET and sum of all PERTURBAGENS affecting all valid KINASES that are potentially affecting a TARGET

% Average inhibition score of all PERTURBAGENS affecting a KINASE on a TARGET	
k_affects_t_with_avg_score(Kinase, Target, Cell_Line, Score) :-
    findall(Score1, pkt_1(_, Kinase, Target, Cell_Line, Score1), ScoreList), % ScoreList is the list of all inhibition scores
    sum(ScoreList, Sum), % Sum is sum of all inhibition scores
    length(ScoreList, Len), % Len is the length of the ScoreList
    (Len is 0 -> Score is 0 ; Score is Sum / Len). % Average inhibition score of all PERTURBAGENS affecting a KINASE on a TARGET is Score
    
% Sum is the sum of all PERTURBAGENS affecting all valid KINASES that are potentially affecting a TARGET
all_k_that_affect_t_with_sum_of_avg_score([], _, _, _, 0).
all_k_that_affect_t_with_sum_of_avg_score([H|T], Kinase, Target, Cell_Line, Sum) :-
   all_k_that_affect_t_with_sum_of_avg_score(T, Kinase, Target, Cell_Line, Rest),
   k_affects_t_with_avg_score(H, Target, Cell_Line, Single_K_Score),
   (H = Kinase -> Similarity is 0 ; k1_k2_is_similar(H, Kinase, _, _, Similarity)),
   Sum is (Single_K_Score * (1 - Similarity)) + Rest.
   
% Perturbagen > Kinase association   
pk(Perturbagen, Kinase, Sum_of_inhibition, Sum_of_excitation, Score) :-
   distinct(pk_has_interaction(Perturbagen, Kinase)), % One PERTURBAGEN > KINASE > TARGET association
   findall(Score1, pkt_2(Perturbagen, Kinase, _, down, _, _, _, _, Score1), Inhibition_list), % Inhibition_list is a list of the scores of all PERTURBAGEN > KINASE > TARGET inhibitions
   findall(Score2, pkt_2(Perturbagen, Kinase, _, up, _, _, _, _, Score2), Excitation_list), % Excitation_list is a list of the scores of all PERTURBAGEN > KINASE > TARGET excitations
   sum(Inhibition_list, Sum_of_inhibition), % Sum_of_inhibition is the sum of all scores in Inhibition_list
   sum(Excitation_list, Sum_of_excitation), % Sum_of_excitation is the sum of all scores in Excitation_list
   (Sum_of_inhibition < 0.0000001 -> Score is 0 ; Score is Sum_of_inhibition * (Sum_of_inhibition / (Sum_of_inhibition + Sum_of_excitation))).

% A PERTURBAGEN and a KINASE has an interaction (used in pk/5)
pk_has_interaction(Perturbagen, Kinase) :- 
   expressedin(Kinase, _), % KINASE iS expressed in at least one cell line
   perturbs(Perturbagen, Target, down, _), % A PERTURBAGEN inhibits a TARGET with a score of
   knowntarget(Target, Kinase). % TARGET is associated with (target of) KINASE

%-----------------------------------------------------------------------------------------------------------------------------------%
% Two KINASES are similar if they have common PERTURBAGENS affecting them
:- table k1_k2_is_similar_by_perturbagens/5.
k1_k2_is_similar_by_perturbagens(K1, K2, Shared_count, Total_unique_perturbagen_count, Ratio) :-
	% kinase(K1),
	% kinase(K2), 
	K1 \= K2,	
	(knowninhibitor(_, K1, _, _) -> findall(Perturbagen1, knowninhibitor(Perturbagen1, K1, _, _), PerturbagenList1) ; PerturbagenList1 = []), % PerturbagenList1 is the list of all PERTURBAGENS that are known to inhibit KINASE1
	list_to_set(PerturbagenList1, PerturbagenList1_unique),	% Remove duplicates
	(knowninhibitor(_, K2, _, _) -> findall(Perturbagen2, knowninhibitor(Perturbagen2, K2, _, _), PerturbagenList2) ; PerturbagenList2 = []), % PerturbagenList2 is the list of all PERTURBAGENS that are known to inhibit KINASE2
	list_to_set(PerturbagenList2, PerturbagenList2_unique), % Remove duplicates
	count2(PerturbagenList1_unique, PerturbagenList2_unique, Shared_count), % Shared_count is the number of common PERTURBAGENS between PerturbagenList1 and PerturbagenList2
	append(PerturbagenList1_unique, PerturbagenList2_unique, Total_perturbagenList), % Add both lists together
	list_to_set(Total_perturbagenList, Total_unique_perturbagenList), % Remove duplicates
	length(Total_unique_perturbagenList, Total_unique_perturbagen_count), % Total_unique_perturbagen_count is the number of PERTURBAGENS that are affecting either KINASE1 or KINASE2
	((Shared_count is 0 ; Total_unique_perturbagen_count = 0) -> Ratio is 0 ; Ratio is Shared_count / Total_unique_perturbagen_count). % Ratio is the ratio of common PERTURBAGENS and the number of PERTURBAGENS that are affecting either KINASE1 or KINASE2

% Two KINASES are similar if they have common TARGETS that they affect
:- table k1_k2_is_similar_by_targets/5.
k1_k2_is_similar_by_targets(K1, K2, Shared_count, Total_unique_target_count, Ratio) :-
	% kinase(K1),
	% kinase(K2),
	K1 \= K2,
	(knowntarget(_, K1) -> findall(Target1, knowntarget(Target1, K1), TargetList1) ; TargetList1 = []), % TargetList1 is the list of all TARGETS that are known to be affected by KINASE1
	list_to_set(TargetList1, TargetList1_unique), % Remove duplicates
	(knowntarget(_, K2) -> findall(Target2, knowntarget(Target2, K2), TargetList2) ; TargetList2 = []), % TargetList2 is the list of all TARGETS that are known to be affected by KINASE2
	list_to_set(TargetList2, TargetList2_unique), % Remove duplicates
	count2(TargetList1_unique, TargetList2_unique, Shared_count), % Shared_count is the number of common TARGETS between TargetList1 and TargetList2
	append(TargetList1_unique, TargetList2_unique, Total_targetList), % Add both lists together
	list_to_set(Total_targetList, Total_unique_targetList), % Remove duplicates
	length(Total_unique_targetList, Total_unique_target_count), % Total_unique_target_count is the number of TARGETS that are affected by either KINASE1 or KINASE2
	((Shared_count is 0 ; Total_unique_target_count = 0) -> Ratio is 0 ; Ratio is Shared_count / Total_unique_target_count). % Ratio is the ratio of common TARGETS and the number of TARGETS that are affected by either KINASE1 or KINASE2

% Two KINASES are similar if they have both common PERTURBAGENS affecting them and common TARGETS that they affect
:- table k1_k2_is_similar/5.
k1_k2_is_similar(K1, K2, P_Score, T_Score, Score) :-
	% kinase(K1),
	% kinase(K2),
	k1_k2_is_similar_by_perturbagens(K1, K2, _, _, P_Score), % Two KINASES are similar if they have common PERTURBAGENS affecting them 
	k1_k2_is_similar_by_targets(K1, K2, _, _, T_Score), % Two KINASES are similar if they have common TARGETS that they affect 
	Score is P_Score * T_Score.
	
%-----------------------------------------------------------------------------------------------------------------------------------%
% Util
printlist([]).
printlist([H|T]) :-
   write(H), write('\n'),
   printlist(T).
   
unique(Kinase, Target, Score, Single_K_Score) :-
   pkt_1(_, Kinase, Target, _, _),
   kt_assoc(Kinase, Target, _, Score, Single_K_Score).

kinase_affects_target(Kinase, Target, Score):-
	knowntarget(Target, Kinase),
	k_affects_t_with_avg_score(Kinase, Target, _, Score).
	
sum([], 0).
sum([H|T], Sum) :-
   sum(T, Rest),
   Sum is H + Rest.
   
count2([],[],0).
count2([],[_|_],0).
count2([H|T],List,S):-
    count(H, List, N),
    count2(T, List, M),
    S is N + M.

count(_, [], 0).
count(X, [X | T], N) :-
  !, 
  count(X, T, N1),
  N is N1 + 1.
  
count(X, [_ | T], N) :-
  count(X, T, N).
	
%-----------------------------------------------------------------------------------------------------------------------------------%
pathway(List):-
	retractall(seen(_)),
	go_up_pathway('AKT1', 'S473', List).


go_up_pathway(Lvl3_kinase, Lvl3_phosphosite, List) :-
	assert(seen(Lvl3_kinase)),
	kt_related(Lvl2_kinase, Lvl3_kinase, Lvl3_phosphosite ,_ ,_, _),
	kt_related(Lvl1_kinase, Lvl2_kinase, Lvl2_phosphosite,_ ,_, _),
	Lvl1_kinase \= Lvl2_kinase,
	Kinase_and_Phospho = [Lvl2_kinase, Lvl2_phosphosite],
	append(List, Kinase_and_Phospho, List2),
	not(seen(Lvl2_kinase)),
	go_up_pathway(Lvl2_kinase, Lvl2_phosphosite, List2).








