function xn=thirty_two2two(yn)
bit_number = length(yn)*5;
xn = zeros(bit_number,1);
bound = 0.5:1:30.5;
for i=1:1:length(yn)
    if(yn(i)<bound(1))
        xn(5*i-4:5*i)=[0,0,0,0,0];
    elseif(yn(i)<bound(2))
        xn(5*i-4:5*i)=[0,0,0,0,1];
    elseif(yn(i)<bound(3))
        xn(5*i-4:5*i)=[0,0,0,1,0];
    elseif(yn(i)<bound(4))
        xn(5*i-4:5*i)=[0,0,0,1,1];
    elseif(yn(i)<bound(5))
        xn(5*i-4:5*i)=[0,0,1,0,0];
    elseif(yn(i)<bound(6))
        xn(5*i-4:5*i)=[0,0,1,0,1];
    elseif(yn(i)<bound(7))
        xn(5*i-4:5*i)=[0,0,1,1,0];  
    elseif(yn(i)<bound(8))
        xn(5*i-4:5*i)=[0,0,1,1,1];
    elseif(yn(i)<bound(9))
        xn(5*i-4:5*i)=[0,1,0,0,0];
    elseif(yn(i)<bound(10))
       xn(5*i-4:5*i)=[0,1,0,0,1];
    elseif(yn(i)<bound(11))
        xn(5*i-4:5*i)=[0,1,0,1,0];
    elseif(yn(i)<bound(12))
        xn(5*i-4:5*i)=[0,1,0,1,1];
    elseif(yn(i)<bound(13))
        xn(5*i-4:5*i)=[0,1,1,0,0];
    elseif(yn(i)<bound(14))
        xn(5*i-4:5*i)=[0,1,1,0,1];
    elseif(yn(i)<bound(15))
        xn(5*i-4:5*i)=[0,1,1,1,0];
    elseif(yn(i)<bound(16))
        xn(5*i-4:5*i)=[0,1,1,1,1];
    elseif(yn(i)<bound(17))
        xn(5*i-4:5*i)=[1,0,0,0,0];
    elseif(yn(i)<bound(18))
        xn(5*i-4:5*i)=[1,0,0,0,1];
    elseif(yn(i)<bound(19))
        xn(5*i-4:5*i)=[1,0,0,1,0];
    elseif(yn(i)<bound(20))
        xn(5*i-4:5*i)=[1,0,0,1,1];  
    elseif(yn(i)<bound(21))
        xn(5*i-4:5*i)=[1,0,1,0,0];
    elseif(yn(i)<bound(22))
        xn(5*i-4:5*i)=[1,0,1,0,1];
    elseif(yn(i)<bound(23))
        xn(5*i-4:5*i)=[1,0,1,1,0];
    elseif(yn(i)<bound(24))
        xn(5*i-4:5*i)=[1,0,1,1,1];
    elseif(yn(i)<bound(25))
       xn(5*i-4:5*i)=[1,1,0,0,0];
    elseif(yn(i)<bound(26))
        xn(5*i-4:5*i)=[1,1,0,0,1];
    elseif(yn(i)<bound(27))
        xn(5*i-4:5*i)=[1,1,0,1,0];
    elseif(yn(i)<bound(28))
        xn(5*i-4:5*i)=[1,1,0,1,1];
    elseif(yn(i)<bound(29))
        xn(5*i-4:5*i)=[1,1,1,0,0];
    elseif(yn(i)<bound(30))
        xn(5*i-4:5*i)=[1,1,1,0,1];
    elseif(yn(i)<bound(31))
        xn(5*i-4:5*i)=[1,1,1,1,0];
    else
        xn(5*i-4:5*i)=[1,1,1,1,1];
    end
end

