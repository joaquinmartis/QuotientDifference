Program QuotientDifference

!1- Ingresa el grado del polinomio
!2- Ingresa los coeficientes de forma ordenada
!3- Comienza QD
    !i-Ejecuta primeras iteraciones

IMPLICIT NONE

integer:: Grado, i ,Metodo,maxiter
real(8), allocatable, dimension(:) :: VecCoeficientes,VecRaices
complex(8), allocatable, dimension(:):: VecRaicesCOMPLEJAS
real(8) tol,diferencia

! COMIENZA EL PROGRAMA

CALL LeeCoeficientes(VecCoeficientes,Grado)
!----------ERROR--------!
write(*,*)'Ingrese la cota maxima de error'
read(*,*)tol
!------Iteraciones-----!
write(*,*)"Ingrese maxima cantidad de iteraciones" 
read(*,*)maxiter

Metodo=0
CALL PreProcesamiento(VecCoeficientes,Grado,Metodo,diferencia)
do while (CoeficienteCero(VecCoeficientes,Grado) .AND. Metodo<4) 

	Metodo=Metodo+1
	CALL PreProcesamiento(VecCoeficientes,Grado,Metodo,Diferencia) !hay 3 metodos para corregir 0s, haciendo diferentes reempazos s=x+d s+cx s=1/x
	!ademas puede bajar de grado si el ultimo termino es 0
	
enddo

if (CoeficienteCero(VecCoeficientes,Grado)) then !si el metodo es 4 enotnces no fue posible alterar el polinomio

  write(*,*)"Imposible obtener raices del polinomio : No fue posible modificarlo para que sus coeficientes fueran distintos de 0"
  
 else
	select case (grado)
	case (1) !grado=0
		write(*,*)"Error grado=0"
	case (2)
		write(*,*)"Raiz: ",VecCoeficientes(1)/VecCoeficientes(2)
	case (3)
		!call Resolvente(VecCoeficientes)!Falta desarrollar
	case default
		call QuotientDifferenceSub(VecCoeficientes,VecRaices,VecRaicesCOMPLEJAS,grado,tol,maxiter)
		!call CorrigeSolucion(VecRaices,grado,Metodo,diferencia)
	end select
endif




CONTAINS

!----------------------------------!
subroutine PreProcesamiento(VecCoeficientes,Grado,Metodo,diferencia)
real(8),intent(inout) :: VecCoeficientes(:)
real(8), intent(out) :: diferencia
integer, intent(inout) ::Metodo,Grado

diferencia = 1 ! variable para el metodo 1 x-d  darle algun valor 
select case(Metodo)
	case (0) !Caso 0 en ultimo lugar
		do while (VecCoeficientes(grado)==0)
			grado=grado-1
		enddo
	case (1)
		call Metodo1(VecCoeficientes,diferencia,Grado)
	case (2)
	case (3)
end select

endsubroutine
!----------------------------------!
subroutine Metodo1(VecCoeficientes,d,Grado)
real(8) ,intent(inout) :: VecCoeficientes(:),d
integer, intent(in) :: Grado
integer i,j,k

do i=2,grado

 k=i
 if(VecCoeficientes(i)/=0) then
   do j=1,i-1
    VecCoeficientes(j) = VecCoeficientes(j) + Combinatoria(i-1,k-1)*(VecCoeficientes(i)**(i-k))*(-1*d)**(k-1)
    k=k-1
   end do
 
  VecCoeficientes(j+1) = (VecCoeficientes(i)**(i-1))
 end if

end do

endsubroutine Metodo1
!----------------------------------!
function Combinatoria(m,n)
integer, intent(in) :: m,n
real(8) Combinatoria

Combinatoria=Factorial(m)/((Factorial(n)*Factorial(m-n))*1.0) !el resultado de la operacion es un entero asi q *1.0 para q sea real

endfunction Combinatoria
!----------------------------------!
function Factorial(x)
integer, intent(in) :: x
integer Factorial

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
real(8) VecCoeficientes(:)
integer Grado,i

do while ((VecCoeficientes(i)/=0) .AND. (i<Grado))
	i=i+1
enddo

CoeficienteCero=((Grado==i).AND.(VecCoeficientes(grado)/=0))

endfunction

!----------------------------------!
SUBROUTINE LeeCoeficientes(VecCoeficientes,Grado)
integer i,Grado
real(8), allocatable :: VecCoeficientes(:)

open(2,FILE="CoeficientesPolinomio.txt", action="read")
read(2,*)grado
grado=grado+1
allocate(VecCoeficientes(grado))

do i=1,grado 
	read(2,*)VecCoeficientes(i)
enddo



END SUBROUTINE
!----------------------------------!

subroutine QuotientDifferenceSub(VecCoeficientes,VecRaices,VecRaicesCOMPLEJAS,grado,tol,maxiter)
!parametros 
real(8) ,intent(inout)::VecCoeficientes(:)
real(8) ,allocatable, intent(inout) :: VecRaices(:)
complex(8), dimension(1:grado), intent(inout) :: VecRaicesCOMPLEJAS
integer ,intent(in):: grado,maxiter
real(8) ,intent(in)::tol
!var locales
integer i,iter
real(8) error
real(8), dimension(1:grado) :: VecRaices_ANT
real(8),dimension(0:grado) :: VecError
real(8) CoefV,CoefU

complex(8) RaizComp1,RaizComp2
!VECError va de 0 a grado
!VEcRaiz va de 1 a grado
!VecCoeficientes va de 1 a grado

   
!COMIENZA QD
CALL CondicionesIniciales(VecCoeficientes,VecRaices,VecError,grado) 
allocate(VecRaices(grado))
iter=1
error=2*tol 

do while (error>tol .AND. maxiter>iter)

	VecRaices_ANT=VecRaices
    do i=1, grado    !Calculo la fila de entera en cada iteración y el grado +1 porque el vector comiena en 1 , el el coef de maximo grado esta een grado +1
       VecRaices(i)=VecError(i)-VecError(i-1)+VecRaices(i)     
    enddo
    
    Do i=1, grado-1      !Como e(0)=e(n)=0 siempre, uso esos límites en el ciclo
       VecError(i)=VecError(i)*(VecRaices(i+1)/VecRaices(i))   ! Hay un error de tipeo en el pdf que figura q(i-1) pero en la tabla utiliza q(i+1), sino no da
    end do
    

     iter=iter+1
     error=maxval(abs(VecError))
end do 
    if (maxiter<iter) then !si salio por iteraciones del do while significa que no se aproxima a 0 el error maximo
		do i=0,grado-1
			if (abs(VecError(i))>tol) then
				CoefU=VecRaices(i)+VecRaices(i+1)
				CoefV=(-1)*VecRaices_ANT(i)*VecRaices(i+1) !No se si va el -1 pero lo vi asi
				CALL MetodoBairstow (VecCoeficientes,grado,tol,CoefU,CoefV)
				CALL Resolvente(CoefU, CoefV, RaizComp1, RaizComp2)
				VecRaicesCOMPLEJAS(i)=RaizComp1
				VecRaicesCOMPLEJAS(i+1)=RaizComp2 ! en i+1
			endif
		enddo
	endif



end subroutine

!----------------------------------!

subroutine CondicionesIniciales(VecCoeficientes,VecRaices,VecError,grado)
real(8), intent(inout)::VecCoeficientes(:),VecRaices(:)
integer ,intent(in):: grado
real(8), intent(inout) :: VecError(0:)

VecRaices=0
VecRaices(1)=-1*(VecCoeficientes(grado-1))/(VecCoeficientes(grado))

do i=1, grado

      VecError(i)=VecCoeficientes(grado-i-1)/VecCoeficientes(grado-i)
      
end do

VecError(0)=0
VecError(grado)=0

endsubroutine
!----------------------------------!

subroutine MetodoBairstow (VecCoeficientes,grado,tol,CoefU,CoefV)
real(8), intent(in) :: VecCoeficientes(:)
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
	do i=grado, 2, -1
	
		Q = VecCoeficientes(i) + CoefU*Q_ant1 + CoefV*Q_ant2
		Q_ant2 = Q_ant1
		Q_ant1 = Q
		P = Q + CoefU*P_ant1 + CoefV*P_ant2
		P_ant3 = P_ant2
        P_ant2 = P_ant1
        
	enddo
	q = VecCoeficientes(1) + CoefU* Q_ant1 + CoefV* Q_ant2

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

end subroutine


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
