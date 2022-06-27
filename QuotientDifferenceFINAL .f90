Program QuotientDifference

!1- Ingresa el grado del polinomio
!2- Ingresa los coeficientes de forma ordenada
!3- Comienza QD
    !i-Ejecuta primeras iteraciones

IMPLICIT NONE
!----------------------
integer ,parameter:: grado=6
!---------------------------
integer maxiter
real(8) tol
real(8), dimension(0:grado) :: VecCoeficientes 
complex(8), dimension(1:grado) :: VecRaices

character res

!-------------
CALL LeeCoeficientes(VecCoeficientes,Grado)
!----------ERROR--------!
!write(*,*)'Ingrese la cota maxima de error'
!read(*,*)tol
tol=0.00001
!------Iteraciones-----!
!write(*,*)"Ingrese maxima cantidad de iteraciones" 
!read(*,*)maxiter
maxiter=30

if (CoeficienteCero(VecCoeficientes,Grado) .eqv. .FALSE.) then
	write(*,*)"Desea realizar algun tipo de tranformacion S/N"
		read(*,*)res
	if ((res=='S' .OR. res=='s') .OR. (CoeficienteCero(VecCoeficientes,Grado))) then
		call MenuPreProcesamiento(VecCoeficientes,Grado,tol,maxiter)
	else
		call QuotientDifferenceSub(VecCoeficientes,VecRaices,grado,tol,maxiter)	
		call MuestraRaices(VecRaices,grado)
	endif
else
	write(*,*)"El polinomio ingresado posee 0 en algun coeficiente, se procedera con una transformacion del mismo"
	call MenuPreProcesamiento(VecCoeficientes,Grado,tol,maxiter)
endif



CONTAINS
!----------------------------------!
subroutine MenuPreProcesamiento(VecCoeficientes,Grado,tol,maxiter)
real(8), intent(inout) :: VecCoeficientes(:)
integer, intent(in) :: Grado,maxiter
real(8), intent(in) :: tol



real(8) constante 
integer Metodo
write(*,*)"Menu de Preprocesamiento del polinomio"
if (CoeficienteCero(VecCoeficientes,Grado)) then
	write(*,*)"El unico metodo disponible para 0s en algun coeficiente es la sutitucion s=(x-d)"
		Metodo=1
else
	write(*,*)"Ingrese tipo de preprocesamiento de 0s"
	write(*,*)"1. s=(x-d)"
	write(*,*)"2. s=cx" 
	write(*,*)"3. s=1/x"
		read(*,*)Metodo
endif


select case(Metodo)
	case (1)
		write(*,*)"El metodo elegido es s=x-d) ingrese valor para d: "
			read(*,*)constante
		call Metodo_xmenosd(VecCoeficientes,constante,Grado)
		
		write(*,*)"Nuevo Polinomio"
		call MuestraPolinomio(VecCoeficientes,grado)
		
		do while (CoeficienteCero(VecCoeficientes,Grado)) !si sigue habiendo 0s elige nuevo valor para d
			constante=-constante
			call Metodo_xmenosd(VecCoeficientes,constante,Grado) ! -diferencia para volver al polinomio original
			
			write(*,*)"El polinomio todavia posee 0s, ingrese un nuevo valor para d"
				read(*,*)constante
				
			call Metodo_xmenosd(VecCoeficientes,constante,Grado)	
			
			write(*,*)"Nuevo Polinomio"
			call MuestraPolinomio(VecCoeficientes,grado)
		enddo
		
		call QuotientDifferenceSub(VecCoeficientes,VecRaices,grado,tol,maxiter)	
		call CorrijeRaices_xmenosd(VecRaices,grado,constante)
		call MuestraRaices(VecRaices,grado)
	case (2)
		write(*,*) "El metodo elegido es s=cx) ingrese valor para c: "
			read(*,*)constante
		call Metodo_cx(VecCoeficientes,grado,constante) 
		call QuotientDifferenceSub(VecCoeficientes,VecRaices,grado,tol,maxiter)
		call CorrigeRaices_cx(VecRaices,grado,constante)
		call MuestraRaices(VecRaices,grado)
	case (3)
		write(*,*)"El metodo elegido es s=1/x "
		call Metodo_1divx(VecCoeficientes,Grado)
		call QuotientDifferenceSub(VecCoeficientes,VecRaices,grado,tol,maxiter)
		call CorrigeRaices_1divx(VecRaices,grado)
		call MuestraRaices(VecRaices,grado)
end select


endsubroutine
!----------------------------------!
subroutine CorrijeRaices_xmenosd(VecRaices,grado,constante)
complex(8) ,intent(inout) :: VecRaices(:)
integer ,intent(in) :: grado
real(8) ,intent(in) :: constante

integer i

do i=1,grado
	VecRaices(i)=VecRaices(i)-constante
enddo

end subroutine

!----------------------------------!

subroutine CorrigeRaices_cx(VecRaices,grado,constante)
complex(8) ,intent(inout) :: VecRaices(:)
integer ,intent(in) :: grado
real(8) ,intent(in) :: constante

integer i

do i=1,grado
	VecRaices(i)=VecRaices(i)*(constante*1.0)
enddo

end subroutine


!----------------------------------!

subroutine CorrigeRaices_1divx(VecRaices,grado)
complex(8) ,intent(inout) :: VecRaices(:)
integer ,intent(in) :: grado

integer i

do i=1,grado
	VecRaices(i)=1./VecRaices(i)
enddo

end subroutine

!----------------------------------!
subroutine MuestraPolinomio(VecPolinomio,grado)
real(8), intent(in) :: VecPolinomio(0:grado)
integer , intent(in) :: grado

integer i


	write(*,'(F7.2)',ADVANCE='NO')VecPolinomio(0)


do i=1,grado

	
		WRITE(*,'(A,F12.2,A,I2)',ADVANCE='NO')" + ",VecPolinomio(i),"x^",i
	
	
enddo
write(*,*)
endsubroutine

!----------------------------------!


subroutine MuestraRaices(VecRaices,grado)
complex(8), intent(in) :: VecRaices(:)
integer i,grado

WRITE(*,*)
WRITE(*,*)

do i=1, grado
		if (imag(VecRaices(i))==0) then
		     write(*,'(A,I2,A,F7.2,A)')"Raiz ",i,": ",real(VecRaices(i))
		else 
			write(*,'(A,I2,A,F7.2,A,F7.2,A)')"Raiz ",i,": ",real(VecRaices(i))," + ",imag(VecRaices(i)),"i "
		endif
enddo



end subroutine

!----------------------------------!

subroutine Metodo_xmenosd(coeficientes,d,grado)
real(8) ,intent(inout) :: coeficientes(0:grado),d
real(8) :: VecCoeficientes(grado+1)
integer, intent(in) :: grado
integer i,j,k

do i=0,grado
VecCoeficientes(i+1) = coeficientes(i)
end do

do i=2,grado+1

 k=i
 if(VecCoeficientes(i)/=0) then
   do j=1,i-1
    VecCoeficientes(j) = VecCoeficientes(j) + VecCoeficientes(i)*Combinatoria(i-1,k-1)*(-1*d)**(k-1.)
    k=k-1
   end do
 
 end if

end do



do i=0,grado
 coeficientes(i) = VecCoeficientes(i+1) 
end do

endsubroutine 


subroutine Metodo_xmenosdBIS(coeficientes,d,grado)
real(8) ,intent(inout) :: coeficientes(0:grado),d
real(8) :: VecCoeficientes(1:grado+1)
integer, intent(in) :: grado
integer i,j,k

do i=0,grado
VecCoeficientes(i+1) = coeficientes(i)
end do

do i=1,grado

 k=i
 if(VecCoeficientes(i)/=0) then
   do j=1,i-1
    VecCoeficientes(j) = VecCoeficientes(j) + VecCoeficientes(i)* Combinatoria(i-1,k-1) * (-1*d)**(k-1.)
    k=k-1
   end do
 
 end if

end do


write(*,*)coeficientes
do i=0,grado
 coeficientes(i) = VecCoeficientes(i+1) 
end do

endsubroutine 
!----------------------------------!
subroutine Metodo_1divx (VecCoeficientes,grado) 
real(8), intent(inout) :: VecCoeficientes(0:)
integer, intent(in) :: grado
 
integer i
real(8) aux
write(*,*) VecCoeficientes

do i=0,grado/2  
	aux=VecCoeficientes(i)
	VecCoeficientes(i)=VecCoeficientes(grado-i)
	VecCoeficientes(grado-i)=aux
enddo

write(*,*) VecCoeficientes

endsubroutine

!----------------------------------!

subroutine Metodo_cx(VecCoeficientes,grado,constante) 
real(8), intent(inout) :: VecCoeficientes(0:)
integer, intent(in) :: grado
real(8), intent(in) :: constante

integer i

do i=0,grado
	VecCoeficientes(i)=VecCoeficientes(i)*(constante**i)
enddo

end subroutine
!----------------------------------!
function Combinatoria(m,n)
integer, intent(in) :: m,n
real(8) Combinatoria

Combinatoria=Factorial(m)/((Factorial(n)*Factorial(m-n))*1.0)

endfunction Combinatoria
!----------------------------------!
function Factorial(x)
integer, intent(in) :: x
integer Factorial
integer i 

if (x==0) then
	Factorial=1
else
	Factorial=1
	do i=1,x
		Factorial=Factorial*i
	enddo
end if	

end function Factorial
!----------------------------------!

function CoeficienteCero(VecCoeficientes,Grado)
logical CoeficienteCero
real(8) VecCoeficientes(0:)
integer Grado,i
i=0
do while ((VecCoeficientes(i)/=0.) .AND. (i<Grado))
	i=i+1
enddo

CoeficienteCero=(VecCoeficientes(i)==0.)

endfunction

!----------------------------------!
SUBROUTINE LeeCoeficientes(VecCoeficientes,Grado)
integer i,Grado
real(8) :: VecCoeficientes(0:)

open(2,FILE="CoeficientesPolinomio.txt", action="read")


write(*,'(A,I2,A)',advance='NO')"Polinomio de grado ",Grado," : "
do i=0,Grado ,1
	read(2,*)VecCoeficientes(i)
	
	!WRITE(*,'(F7.2,A,I2,A)',ADVANCE='NO')VecCoeficientes(i),"x^",i," + "
	
enddo
call MuestraPolinomio(VecCoeficientes,grado)
write(*,*)

END SUBROUTINE
!----------------------------------!

subroutine QuotientDifferenceSub(VecCoeficientes,VecRaices,grado,tol,maxiter)
!parametros 
real(8) ,intent(inout)::VecCoeficientes(:)
complex(8), dimension(1:grado), intent(inout) :: VecRaices
integer ,intent(in):: grado,maxiter
real(8) ,intent(in)::tol
!var locales
real(8), dimension(1:grado) :: VecQ_ANT
real(8) ,dimension(1:grado) :: VecQ
integer i,iter
real(8) error
real(8),dimension(0:grado) :: VecError
real(8) CoefV,CoefU
complex(8) RaizComp1,RaizComp2

   
CALL CondicionesIniciales(VecCoeficientes,VecQ,VecError,grado) 

iter=1
error=2*tol 

write(*,*)"iter",iter
write(*,'(A)',ADVANCE='no')"Q :           "
call MuestraVector(VecQ,grado)
write(*,'(A)',ADVANCE='no')"Errores : "
call MuestraVector(VecError,grado+1)

do while (error>tol .AND. maxiter>iter)
	write(*,*)"-------------------------"
	VecQ_ANT=VecQ
    do i=1, grado                                    !Calculo la fila de entera en cada iteración y el grado +1 porque el vector comiena en 1 , el el coef de maximo grado esta een grado +1
       VecQ(i)=VecError(i)-VecError(i-1)+VecQ(i)     
       
    enddo
    
    
    Do i=1, grado-1                                   !Como e(0)=e(n)=0 siempre, uso esos límites en el ciclo
       VecError(i)=VecError(i)*(VecQ(i+1)/VecQ(i))    ! Hay un error de tipeo en el pdf que figura q(i-1) pero en la tabla utiliza q(i+1), sino no da
    end do
    

     iter=iter+1
     error=maxval(abs(VecError))
     write(*,*)"iter",iter
     write(*,'(A)',ADVANCE='no')"Q :           "
     call MuestraVector(VecQ,grado)
     write(*,'(A)',ADVANCE='no')"Errores : "
     call MuestraVector(VecError,grado+1)
end do 
		do i=0,grado-1
			if (abs(VecError(i))>tol) then
				CoefU=VecQ(i)+VecQ(i+1)
				CoefV=(-1)*VecQ_ANT(i)*VecQ(i+1) !No se si va el -1 pero lo vi asi
				CALL MetodoBairstow (VecCoeficientes,grado,tol,CoefU,CoefV)
				CoefV=-CoefV
				CoefU=-CoefU
				CALL Resolvente(CoefU, CoefV, RaizComp1, RaizComp2)
				VecRaices(i)=RaizComp1
				VecRaices(i+1)=RaizComp2
			else
				VecRaices(i+1)=DCMPLX(VecQ(i+1), 0.)
			endif
		enddo

end subroutine

!----------------------------------!

subroutine CondicionesIniciales(VecCoeficientes,VecRaices,VecError,grado)
real(8), intent(inout)::VecCoeficientes(0:),VecRaices(1:)
integer ,intent(in):: grado
real(8), intent(inout) :: VecError(0:) !NO SACAR ESTE 0, SI LO SACAS REVIENTA EL PROGRAMA
integer i

VecRaices=0.
VecRaices(1)=-1.*(VecCoeficientes(grado-1))/(VecCoeficientes(grado))

do i=1, grado-1

      VecError(i)=VecCoeficientes(grado-i-1)/VecCoeficientes(grado-i)
      
end do

VecError(0)=0.
VecError(grado)=0.

endsubroutine

!----------------------------------!

subroutine MuestraVector(Vector,N)
real(8), intent(in) :: Vector(:)
integer ,intent(in) :: N
integer i


do i=1,N
	write(*,'(F12.3)',ADVANCE='NO')Vector(i)
enddo
write(*,*) 


end subroutine

!----------------------------------!

subroutine MetodoBairstow (VecCoeficientes,grado,tol,CoefU,CoefV)
!parametros
real(8), intent(in) :: VecCoeficientes(0:)
real(8), intent(inout) :: CoefU,CoefV
real(8), intent(in) :: tol
integer, intent(in) ::grado
!var locales
real(8) Q, Q_ant1, Q_ant2 ,P, P_ant1, P_ant2, P_ant3, IncrementoH, IncrementoK, error
integer i

error=2*tol
do while (error >= tol)
	Q_ant1 = 0.
    Q_ant2 = 0.
    P_ant1 = 0.
    P_ant2 = 0.
	do i=grado, 1, -1
	
		Q = VecCoeficientes(i) + CoefU*Q_ant1 + CoefV*Q_ant2
		Q_ant2 = Q_ant1
		Q_ant1 = Q
		
		P = Q + CoefU*P_ant1 + CoefV*P_ant2
		P_ant3 = P_ant2
        P_ant2 = P_ant1
        P_ant1= P
	enddo
	q = VecCoeficientes(0) + CoefU* Q_ant1 + CoefV* Q_ant2
   
	IncrementoH = (Q * P_ant3 - Q_ant1 * P_ant2) / (P_ant2**2. - P_ant1 * P_ant3)
	IncrementoK = (Q_ant1 * P_ant1 - Q * P_ant2) / (P_ant2**2. - P_ant1 * P_ant3)
	
	CoefU = CoefU +IncrementoH
    CoefV = CoefV +IncrementoK

	if(abs(Q) > abs(Q_ant1))then
		error = abs(Q)
	else
		error= abs(Q_ant1)
	endif
	
enddo

endsubroutine

!----------------------------------!

SUBROUTINE Resolvente(b, c, c1, c2)

REAL(8) a
REAL(8), INTENT(IN) :: b, c
COMPLEX(8) Discriminante 
COMPLEX(8), INTENT(OUT) :: c1, c2

    a = 1.
    Discriminante = b**2 - 4. * a * c
    c1 = (-b + sqrt(Discriminante)) / (2. * a)
    c2 = (-b - sqrt(Discriminante)) / (2. * a)

END SUBROUTINE Resolvente


!----------------------------------!
End program
