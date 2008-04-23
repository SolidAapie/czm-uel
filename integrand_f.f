c      xi = z(1)
c      ze = z(2)
      call user_element_f(z(1),z(2),interface_z,area,
     &     thickness_coord,nope,
     &     xl,u_node_interface,u_max_node_interface,properties,
     &     integrand)
      do nfuncs_count=1,nfuncs_
c This is a kludge that uses the integer truncation from (*/3) = floor(*/3)
         f_return(nfuncs_count) = integrand(nfuncs_count-
     &        3*((nfuncs_count+2)/3)+3,(nfuncs_count+2)/3)
      enddo
