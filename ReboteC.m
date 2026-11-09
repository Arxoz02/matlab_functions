function [curva, trayectoria, intersecciones] = reboteC(p0,v0,fx,fy,k,g,R,rt,rs)
    % Calcula la trayectoria de una partícula al rebotar dentro de una curva paramétrica
    % y bajo una influencia gravitacional
    %
    % p0 :: Posición compleja inicial
    % v0 :: Velocidad compleja inicial
    % fx :: Función paramétrica del eje X en función de s
    % fy :: Función paramétrica del eje Y en función de s
    % k  :: Constante elástica (default=1)
    % g  :: Aceleración sobre la partícula (default=-9.81i)
    % R  :: Cantidad de rebotes (default=6 rebotes)
    % rt :: Calidad de la trayectoria (default=1000 nodos)
    % rs :: Calidad de la curva (default=1000 nodos)
    arguments
        p0; v0; fx; fy; k = 1; g = -9.81i; R = 6; rt = 1000; rs = 1000;
    end
    syms s t real
    %-- Preasignar arreglos --%
    T = zeros(1,R);
    P = [p0 T];
    V = [v0 T];

    df = matlabFunction(diff(fy,s)/diff(fx,s)); % Define la derivada

    for n = 1:R
        p = P(n) + V(n)*t + 0.5*g*t.^2; % Trayectoria de la pelota

        [sT, sF] = vpasolve([real(p) == fx; imag(p) == fy], [t; s], [1e-8 inf; -inf inf]); % Intersección en tiempo y espacio.
        if isempty(sT) % Si no encuentra intersección, rompe el ciclo
            n = n-1;
            P = P(1:n);
            V = V(1:n);
            disp(['Sólo se realizaron ' num2str(n) ' rebotes.'])
            break
        end
    
        [~, ind] = min(double(sT)); % Localiza la primera interección temporal
    
        sF = double(sF(ind)); % Localiza la primera interección espacial
    
        T(n) = double(sT(ind)); % Registro de las intersecciones temporales
    
        O = atan(df(sF)); % Ángulo de la pendiente en la intersección
    
        P(n+1) = P(n) + V(n)*T(n) + 0.5*g*T(n)^2; % Siguiente posición inicial
    
        V(n+1) = V(n) + g*T(n); % Siguiente velocidad inicial
    
        V(n+1) = k*exp(O*1i)*conj(V(n+1)*exp(-1i*O)); % Redirección de la velocidad por el rebote
    end

    %-- Creación de la trayectoria
    TT = linspace(0,sum(T),rt); % Arreglo temporal
    trayectoria = zeros(1,rt); % Preasignación de trayectorias
    for j = 1:n % Se elige 'n' y no 'R' por si se genera algún error durante la búsqueda de intersecciones
        if j == 1 % Para la primera intersección
            ti = (0 <= TT) & (TT <= T(j));
            Ti = TT(ti);
        else % Para las demás intersecciones
            ti = (sum(T(1:j-1)) < TT) & (TT <= sum(T(1:j)));
            Ti = TT(ti)-sum(T(1:j-1));
        end
        % ti :: Busca los intervalos de tiempo de las trayectorias
        % Ti :: Arreglo temporal de la j-ésima trayectoria
        trayectoria(ti) = P(j) + V(j)*Ti + 0.5*g*Ti.^2; % Asigna trayectoria
    end
    intersecciones = P;
    
    %-- Creación de la curva
    C = linspace(-10,10, rs);
    curva = double(subs(fx,C)+1i*subs(fy,C));
end
