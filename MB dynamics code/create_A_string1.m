function [b] = create_A_string1(a)
    b = "[" + replace(join(replace(string(a),[" - ","- "," + ",],["-","-","+"]))," ",", ") + "]";
    b = replace(string(b),["-","+"],[" - "," + ",]);
end