

function loadsdfiter(data::String, no_halt::Bool, assign_descriptors::Bool)
    linesplit
    mol_supplier(data, false, assign_descriptors)
end


function loadsdfiter(file::IOStream, no_halt::Bool, assign_descriptors::Bool)
    data = read(file, String)
    decode
    loadsdfiter(data, false, assign_descriptors)
end


function loadsdfmol(data::String, assign_descriptors::Bool)
    moliter = loadsdfiter(data, false, assign_descriptors)
    next(moliter)
end


function loadsdfmol(file::IOStream, assign_descriptors::Bool)
    data = read(file, String)
    decode
    loadsdfmol(data, assign_descriptors)
end


function sdfblock(lines::Array)
    for line in lines
        if startswith(line, "$$$$")
            yield
        end
        
    end
end


function molsupplier(lines, nohalt, precalc)
    for i, mol, opt in sdfblock(lines)
end
