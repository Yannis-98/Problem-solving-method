module parametre_mod
	real,    parameter :: pi=ACOS(-1.), eps=1e-20
end module parametre_mod

module fonction_pb
contains        
        real function phi(xi) ! solution exacte !
        use parametre_mod 
        implicit none
                real, intent(in) :: xi
                real ::  alpha
                alpha=(7-3*exp(2*pi))/(10.0*(exp(-pi)-exp(2*pi)))        
                phi=alpha*(exp(-xi)-exp(2*xi))+(3.0/10)*exp(2*xi)-(3/10.0)*cos(xi)-(1/10.0)*sin(xi)
        end function phi

        real function fi(xi) !fonction f(x) !
	        real, intent(in) :: xi       
                fi=-cos(xi)
        end function fi

end module fonction_pb

module initialisation_mod
contains
        subroutine initialisation(nx,h,a1,a2,a3,w) !initialisation des paramètre de notre equation discrete !
        use parametre_mod
        implicit none
                integer, intent(in)  :: nx
                real,    intent(out) :: h,a1,a2,a3,w 
                h=pi*1.0/nx
                a1=1+h/2.0 
                a2=2+2*h*h
                a3=-1+h/2.0
                w=1.5
        end subroutine initialisation

        subroutine initialisation_maillage(h,x)
        implicit none
                real, intent(in) :: h
                real, dimension(0:), intent(out) :: x
                integer :: nx,i
                nx=size(x)-1
                do i=0, nx !initialisation de x!
                        x(i)=i*h
                end do
        end subroutine initialisation_maillage

	subroutine initialisation_vecteurs(Aw, Ap, Ae, b, h, a1, a2,a3, x) !initialisation des vecteurs de notre probleme!
		use fonction_pb
		implicit none
		integer :: i, nx
		
		real, dimension(:), allocatable 	:: x, b,Aw, Ae, Ap
		REAL :: a1,a2,a3, h
		!Initialisation de nos vecteur !
		nx=size(x)-1
		DO i=1, (nx-1), 1
			Aw(i)=-a1
			Ap(i)=a2
			Ae(i)=a3
			b(i)=h*h*(fi(x(i)))
		END DO
		b(nx-1)=b(nx-1)-a3
		Aw(1)=0
		Ae(nx-1)=0
	end subroutine initialisation_vecteurs

end module initialisation_mod



module entrees_sorties_mod !ecriture des fichiers!
contains
        subroutine ecriture(nom_fichier,v1,v2, n, format1) 
        use fonction_pb
        implicit none
                character(len=*), intent(in) :: nom_fichier, format1
                real, dimension(0:), intent(in) :: v1
                real, dimension(:), intent(in) :: v2
		integer :: i,nx, n
                nx=size(v1)-1
                open(unit=n,file=nom_fichier)
                        do i=1,nx-1
                                write(unit=n,fmt=format1) v1(i),fi(v1(i)),phi(v1(i)),v2(i)
                        enddo
                close(unit=n)
        end subroutine ecriture

	subroutine ecriture2(nom_fichier,v1,v2, n, format1)
        use fonction_pb
        implicit none
                character(len=*), intent(in) :: nom_fichier, format1
                real, dimension(0:), intent(in) :: v1
                real, dimension(0:), intent(in) :: v2
		integer :: i,nx, n
                nx=size(v1)-1
                open(unit=n,file=nom_fichier)
                        do i=0,nx
                                write(unit=n,fmt=format1) v1(i),v2(i)
                        enddo
                close(unit=n)
        end subroutine ecriture2

end module entrees_sorties_mod

module ordre_methode_mod
contains

	subroutine ordre_methode(nx,u, Ec, E_u, h, x)
		use fonction_pb
		implicit none
		integer :: i
		integer, intent(in) :: nx
		real, intent(in) :: h
		real, dimension(0:), intent(in) :: u, x
		real, dimension(:), intent(out) :: Ec, E_u	
		do i=1, (nx-1) !calcul des erreurs de convergence et consistance!
			Ec(i)=phi(x(i)) - u(i)
			E_u(i)=fi(x(i))+(1./(h*h))*(u(i-1)-2*u(i)+u(i+1))+(1./(2*h))*(-u(i+1)+u(i-1))-2*u(i)
		end do
		open (10, file='consistance_convergence_dec_LU.out')  !stockage dans un fichier!
		do i=1, nx-1
			write(10,fmt='(4e15.6)') E_u(i), Ec(i)
		end do
		close (10)
	end subroutine ordre_methode

end module ordre_methode_mod

module sor_mod
contains
	subroutine sor(s2, nx, r, b, Aw, Us, Ap, Ae, s, dU, Us_k1, w) !algo sor!
		use parametre_mod
		use initialisation_mod
		implicit none
		integer :: i, j
		integer, intent(in) :: nx
		real, intent(inout) :: s2, s
		real, intent(in) :: w
		real, dimension(:), intent(in) :: Ap, Aw, Ae, b
		real, dimension(:), intent(out) :: r
		real, dimension(0:), intent(out) :: Us_k1
		real, dimension(0:), intent(inout) :: Us, dU
		s=0
		s2=1e+20 !parametre permettant d'entrer dans la boucle une premiere fois !
		do i=0, nx ! parametrage des vecteurs !
			Us(i)=0
		        dU(i)=0
			Us_k1(i)=0
		end do
		do i=1, nx-1
			r(i)=0
		end do

		do while ((s2>eps))
			Us(0)=0
			Us(nx)=1
			do j=1, nx-1,  1
				r(j)=b(j)-Aw(j)*Us(j-1)-Ap(j)*Us(j)-Ae(j)*Us(j+1)
				s=(r(j)*r(j))+s
			end do
			s2=sqrt(s)
			if (s2<=eps) then
			else
			do i=1,nx-1
				dU(i)=(w*1.)/(Ap(i))*(r(i)-Aw(i)*dU(i-1))
				Us_k1(i)=Us(i)+dU(i)
			end do
			do i=0, nx, 1
				Us(i)=Us_k1(i)
			end do
			end if
			s=0
		end do
		print *, 'Execution de la methode iterative SOR terminee'
	end subroutine sor

end module sor_mod


PROGRAM convergence
        use parametre_mod
        use fonction_pb
        use initialisation_mod
        use entrees_sorties_mod
	use sor_mod
	use ordre_methode_mod
        implicit none
 
        integer         :: nx=1, i, j, maxi=3 !nb iteration boucle!
	real, dimension(:), allocatable 	:: b, x, y, Aw, Ae, Ap
	real :: s=0,s2=1e+20, w !w est le parametre de relaxation !
	REAL :: h, ti, tf
	REAL :: a1,a2,a3, ninf, const, ninf_conv, const2 !ninf_conv et const2 permettent d'interpreter la convergence !

	REAL, DIMENSION(:), ALLOCATABLE 	:: Lp, Lw, Ue, u, Us, Us_k1, r, dU, Ec, E_u

        do j=1,maxi ! On va jusqu'a nx=10000 car trop long sinon !
		nx=nx*10
		print *, 'on affiche nx ', nx
		!Initialisation des parametres!
		call initialisation(nx,h,a1,a2,a3,w)
		allocate (b(nx-1), x(0:nx), y(nx-1), Aw(nx-1), Ae(nx-1), Ap(nx-1), r(1:nx-1))
		allocate (Lp(nx-1), Lw(nx-1), Ue(nx-1), Us(0:nx), Us_k1(0:nx),dU(0:nx), Ec(nx), E_u(nx))
		!Initialisation de x !
		call initialisation_maillage(h,x)

		!Initialisation de nos vecteurs pour la maillage decale!
		Call initialisation_vecteurs(Aw, Ap, Ae, b, h, a1, a2, a3, x)

		!SOR !
		call cpu_time(ti)
		Call sor(s2, nx, r, b, Aw, Us, Ap, Ae, s, dU, Us_k1, w)
		call cpu_time(tf)
		print *, "Temps d'execution en seconde de la methode SOR : ", tf-ti
		!ecriture de u pour hmax!
		if (j==maxi) then 
			Call ecriture2('Solutions_Usor.out',x,Us, 20, '(4e15.6)')
		end if

		Call ordre_methode(nx,Us, Ec, E_u, h, x) ! calcul des erreurs !

		
		!Interpretation des erreurs de convergence!
		ninf_conv=MAXVAL(ABS(Ec))
		const2=ninf_conv*nx*nx
		print *, 'La norme infini de la convergence vaut :', ninf_conv
		print *, 'La constante de la convergence vaut :', const2
		print *, ''
		open (unit=20, file='convergence_cla_sor.out',position='append')
		        write(unit=20,fmt='(I8, 3e15.6)') nx, h, ninf_conv, const2
		close(unit=20)

		!Interpretation des erreurs de consistance avec u !
		NINF=MAXVAL(ABS(E_u))
		const=NINF*nx*nx
		print *, 'La norme infini de la consistance vaut :', ninf
	    	print *, 'La constante du schema avec u_LU vaut : ', const
		print *, ''
		
		open (unit=20, file='consistance_cla_sor.out',position='append')
		        write(unit=20,fmt='(I8, 3e15.6)') nx, h, ninf, const
		close(unit=20)

		deallocate (b, x, y, Aw, Ae, Ap, r)
		deallocate (Lp, Lw, Ue, Us, Us_k1,dU, Ec, E_u)
		print *, '--------------------'
		print *, ''
	end do
		call execute_command_line('gnuplot -p trace_conv_clas_SOR.gnu') !trace la solution avec le dernier h calculé!
	 
END PROGRAM convergence


	


