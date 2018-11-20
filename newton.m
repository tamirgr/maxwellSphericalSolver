function [root, n] = newton(fn, init, trustr, varargin)

% Newton's method
% with iteration limit
% throws error when converged root is outside specified region

% parameters
tolstop = 1e-14;
tolwarn = 1e-10;
iterlim = 30;

% initialize iteration variables
n = 0;
prev = init;

while true
  root = prev - fn(prev, varargin{:});

  % exit conditions
  diff = abs(root - prev)/abs(prev);
  if diff < tolstop
    break
  elseif n == iterlim;
    break
  else
    prev = root;
    n = n + 1;
  end
end

%n
% errors and warnings
if trustr ~= realmax && (abs(root-init) > trustr || ~isfinite(root))
  error('Trust region violated')
elseif n == iterlim && diff > tolwarn;
  warning('Newton:IterLimitReached', ['Iteration limit was reached. Relative size of final step was ' num2str(diff) '.'])
end
