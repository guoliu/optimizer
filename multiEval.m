function fitness = multiEval(func, inds)
para =  num2cell(inds,2);
fitness = cellfun(func, para);