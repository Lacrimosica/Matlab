function [vector ,n] = extension(int_vector)
    index = 1;
    beginning = 1;
    ending = size(int_vector);
    begin_set=false;
    

    for x = int_vector
        if(x ~= 0.0 && begin_set == true)
            ending = index;
        end

        if(x ~= 0.0)
            if(begin_set == false)
                beginning = index;
                begin_set=true;
            end
           
        end

        index= index +1;
    end
    vector = int_vector(beginning:ending);
    n= ending - beginning +1
end

