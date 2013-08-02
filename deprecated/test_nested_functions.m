
function [addk,getk,gets] = test_nested_functions()
warning('This function is deprecated')


k = 0;

addk = @add_k;
getk = @get_k;
gets = @get_setter;

  function add_k()
  k = k + 1;
  end
  
  function y = get_k()
  y = k;
  end
  
  function f = get_setter(x)
  f = @set_k;
    function set_k()
    k = x;
    end
  end
  
end