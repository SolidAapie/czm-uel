      call user_e_c3d_f(z(1),z(2),interface_z,area,thickness_coord,
     &     nope,xl,u_node_interface,u_max_node_interface,properties,
     &     integrand)
      do nfuncs_count=1,nfuncs_
c This is a kludge that uses the integer truncation from (*/3) = floor(*/3)
         f_return(nfuncs_count) = integrand(nfuncs_count-60*
     &        ((nfuncs_count+59)/60)+60,(nfuncs_count+59)/60)
      enddo
