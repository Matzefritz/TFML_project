%–– what we’re learning
head_pred(belongs,2).

%–– background predicates
body_pred(red,3).
body_pred(gt,2).
body_pred(lt,2).

%–– types
type(belongs,(molecule,atom)).
type(red,(molecule,mult,num)).    % multiplicity is an atom: s, d, t, q
type(gt,(num,num)).
type(lt,(num,num)).

%–– which B-preds may appear in a belongs/2 clause
determination(belongs/2,red/3).
determination(belongs/2,gt/2).
determination(belongs/2,lt/2).

%–– exactly one molecule var per clause
:-
    clause(C),
    #count{V : var_type(C,V,molecule)} != 1.

%–– FOIL’s “one rule, ≤6 literals” bias
max_clauses(1).
max_body(6).
