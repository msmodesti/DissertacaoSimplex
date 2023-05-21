% Load the objective function
fun = @Booth;
esperado=0;

% Define the initial simplex
x0 = [1, 1];  % Initial guess
simplex = [x0; x0 + [1, 0]; x0 + [0, 1]];

% Set the termination criterion
tol = 1e-10;

% Run the Nelder-Mead algorithm
tic;
iter = 0;
while true
    % Sort the simplex based on function values
    values = zeros(size(simplex, 1), 1);
    for i = 1:size(simplex, 1)
        values(i) = fun(simplex(i, :));
    end
    [values, order] = sort(values);
    simplex = simplex(order, :);

    % Calculate the centroid (excluding the worst point)
    centroid = mean(simplex(1:end-1, :));

    % Reflect the worst point through the centroid
    reflected = 2 * centroid - simplex(end, :);

    % Evaluate the reflected point
    reflectedValue = fun(reflected);

    if reflectedValue < values(end-1)
        % Expand
        expanded = 3 * centroid - 2 * simplex(end, :);
        expandedValue = fun(expanded);

        if expandedValue < reflectedValue
            simplex(end, :) = expanded;
            values(end) = expandedValue;
        else
            simplex(end, :) = reflected;
            values(end) = reflectedValue;
        end
    elseif reflectedValue >= values(end-1)
        % Contract
        contracted = 0.5 * (centroid + simplex(end, :));
        contractedValue = fun(contracted);
        if contractedValue < values(end)
            simplex(end, :) = contracted;
            values(end) = contractedValue;
        else
            simplex = 0.5 * (simplex(1, :) + simplex(2:end, :));  % Shrink
            values = zeros(size(simplex, 1), 1);
            for i = 1:size(simplex, 1)
                values(i) = fun(simplex(i, :));
            end
        end
    end

    % Check termination criterion
    if max(abs(values - values(1))) <= tol
        break;
    end

    iter = iter + 1;
end
tempo=toc;

% Display the results
minimizer = simplex(1, :);
minimum = fun(minimizer);
erro = norm(esperado-minimum);
disp('Nelder-Mead Algorithm:');
disp(['Minimum value: ' num2str(minimum)]);
disp(['Minimizer: ' num2str(minimizer)]);
disp(['Iterations: ' num2str(iter)]);
disp(['Erro: ' num2str(erro)]);
disp(['Tempo de máquina:' num2str(tempo)]);

