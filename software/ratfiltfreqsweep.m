% Code to accompany the paper
%
%     Austin, A. P., Horesh, L. and Kalantzis, V.  A rational filtering
%       algorithm for sequences of shifted symmetric linear systems with
%       applications to frequency repsonse analysis.  Submitted, 2023.
%
% Requires Chebfun, the PDE Toolbox, and, optionally, the Parallel Computing
% Toolbox.  Test matrices other than the fem2D and fem3D families may be found
% in the SuiteSparse matrix collection at
%
%     https://sparse.tamu.edu/
%
% Copyright (C) 2023 Anthony P. Austin, Lior Horesh, and Vasileios Kalantzis.
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
%
% 3. Neither the names of the copyright holders nor the names of its contributors
%    may be used to endorse or promote products derived from this software without
%    specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

rng(67714070);
maxNumCompThreads(1);

runTests();
findFEMMeshParameters();
generateFEMMatrices();

%%%%%%%%%% ITERATIONS - FEM PROBLEMS %%%%%%%%%%

problems = {'fem2D-10K',
	    'fem2D-20K',
	    'fem2D-50K',
	    'fem2D-100K',
	    'fem2D-200K',
	    'fem2D-500K',
	    'fem3D-10K',
	    'fem3D-20K',
	    'fem3D-50K',
	    'fem3D-100K',
	    'fem3D-200K',
	    'fem3D-500K'};

opts = struct();
opts.factorization           = 'ldl';
opts.subspace.algorithm      = 'ratfilt';
opts.subspace.tol            = 1.0e-12;
opts.subspace.maxIter        = 50;
opts.subspace.insideEigsOnly = true;
opts.solveForCorrection      = true;
opts.solver.algorithm        = 'gmres';
opts.solver.maxIter          = 50;
opts.solver.restart          = [];
opts.solver.useInitialGuess  = true;
opts.displayLevel            = 2;
opts.parallelize             = true;

matFile = 'iterationExperimentsFEMLowInsideOnly.mat';
outFile = 'iterationExperimentsFEMLowInsideOnly.out';
runIterationExperiments(problems, 'low', opts, matFile, outFile);
printIterationProblemInfoTable(matFile);
computeIterationExperimentsGamma(matFile);
printIterationExperimentsResultsTable(matFile);

opts.subspace.insideEigsOnly = false;

matFile = 'iterationExperimentsFEMLow.mat';
outFile = 'iterationExperimentsFEMLow.out';
runIterationExperiments(problems, 'low', opts, matFile, outFile);
printIterationProblemInfoTable(matFile);
computeIterationExperimentsGamma(matFile);
printIterationExperimentsResultsTable(matFile);

%%%%%%%%%% ITERATIONS - SUITESPARSE PROBLEMS %%%%%%%%%%

problems = {'Kuu/Muu',
	    'crystk03/crystm03',
	    'qa8fk/qa8fm',
	    'windscreen',
	    'bs01'};

opts = struct();
opts.factorization           = 'ldl';
opts.subspace.algorithm      = 'ratfilt';
opts.subspace.tol            = 1.0e-12;
opts.subspace.maxIter        = 50;
opts.subspace.insideEigsOnly = true;
opts.solveForCorrection      = true;
opts.solver.algorithm        = 'gmres';
opts.solver.maxIter          = 50;
opts.solver.restart          = [];
opts.solver.useInitialGuess  = true;
opts.displayLevel            = 2;
opts.parallelize             = true;

matFile = 'iterationExperimentsSSMidInsideOnly.mat';
outFile = 'iterationExperimentsSSMidInsideOnly.out';
runIterationExperiments(problems, 'mid', opts, matFile, outFile);
printIterationProblemInfoTable(matFile);
computeIterationExperimentsGamma(matFile);
printIterationExperimentsResultsTable(matFile);

matFile = 'iterationExperimentsSSLowInsideOnly.mat';
outFile = 'iterationExperimentsSSLowInsideOnly.out';
runIterationExperiments(problems, 'low', opts, matFile, outFile);
printIterationProblemInfoTable(matFile);
computeIterationExperimentsGamma(matFile);
printIterationExperimentsResultsTable(matFile);

opts.subspace.insideEigsOnly = false;

matFile = 'iterationExperimentsSSMid.mat';
outFile = 'iterationExperimentsSSMid.out';
runIterationExperiments(problems, 'mid', opts, matFile, outFile);
printIterationProblemInfoTable(matFile);
computeIterationExperimentsGamma(matFile);
printIterationExperimentsResultsTable(matFile);

matFile = 'iterationExperimentsSSLow.mat';
outFile = 'iterationExperimentsSSLow.out';
runIterationExperiments(problems, 'low', opts, matFile, outFile);
printIterationProblemInfoTable(matFile);
computeIterationExperimentsGamma(matFile);
printIterationExperimentsResultsTable(matFile);

%%%%%%%%%% TIMING - FEM PROBLEM %%%%%%%%%%

opts = struct();
opts.K                       = 16;
opts.factorization           = 'ldl';
opts.subspace.algorithm      = 'ratfilt';
opts.subspace.tol            = 1.0e-12;
opts.subspace.maxIter        = 50;
opts.subspace.insideEigsOnly = false;
opts.solveForCorrection      = true;
opts.solver.algorithm        = 'gmres';
opts.solver.tol              = 1.0e-6;
opts.solver.maxIter          = 50;
opts.solver.restart          = [];
opts.solver.useInitialGuess  = true;
opts.displayLevel            = 2;
opts.parallelize             = false;

matFile = 'timingFEM.mat';
outFile = 'timingFEM.out';
runTimingExperiments('fem3D-500K-T', 'low', opts, matFile, outFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function runIterationExperiments(problems, type, opts, matFile, outFile, resume)
	K = [8, 12, 16, 20, 24];
	tol = [1.0e-8];

	if (nargin < 6)
		resume = false;
	end

	maxNumCompThreads(1);

	diary(outFile);

	if (resume)
		load(matFile);
	else
		data = dictionary();
	end

	for (n = 1:1:numel(problems))
		name = problems{n};
		fprintf('---------- PROBLEM %s (%s) ----------\n', ...
			name, char(datetime('now')));

		[A, M, wmin, wmax, w, f] = loadProblem(name, type);

		if (resume && isKey(data, name))
			dataCell = data(name);
			dataCell = dataCell{1};
		else
			dataCell = cell(numel(K), numel(tol));
		end

		for (i = 1:1:numel(K))
			for (j = 1:1:numel(tol))
				if (resume & ~isempty(dataCell{i, j}))
					fprintf('Found results for K = %d, tol = %.1e\n', K(i), tol(j));
					continue;
				end

				opts.K          = K(i);
				opts.solver.tol = tol(j);
				fprintf('Running K = %d, tol = %.1e\n', ...
				       	K(i), tol(j));

				[u, dataij] = solveShiftedSystems(A, M, f, wmin, wmax, w, opts);

				dataij.opts = opts;
				dataij.type = type;
				dataCell{i, j} = dataij;

				data(name) = {dataCell};
				save(matFile, 'data');
			end
		end
	end

	diary('off');
	save(matFile, 'data');
end

function printIterationProblemInfoTable(matFile)
	load(matFile);

	problems = keys(data);
	maxNameLength = max(cellfun(@length, problems));
	nameFormat = sprintf('%%-%ds', maxNameLength);

	fprintf([nameFormat '       Dim.          a           b    # eigs.\n'], 'Problem');
	fprintf('----------------------------------------------------------------\n');
	for (n = 1:1:numel(problems))
		name = problems{n};
		problemData = data(name);
		problemData = problemData{1};
		[A, W, wmin, wmax, w] = loadProblem(name, problemData{1}.type);

		fprintf([nameFormat '    %6d    % .2e    %.2e        %3d\n'], ...
		        name, size(A, 1), wmin, wmax, ...
		        problemData{1}.subspace.eigCount);
	end
end

function computeIterationExperimentsGamma(matFile)
	load(matFile);

	problems = keys(data);

	opts = struct();
	opts.factorization           = 'ldl';
	opts.subspace.algorithm      = 'ratfilt';
	opts.subspace.tol            = 1.0e-12;
	opts.subspace.maxIter        = 50;
	opts.subspace.insideEigsOnly = false;
	opts.displayLevel            = 2;
	opts.parallelize             = true;

	for (n = 1:1:numel(problems))
		name = problems{n};
		problemData = data(name);
		problemData = problemData{1};

		if (strcmp(name, 'fem3D-500K'))
			opts.K = 16;
		else
			opts.K = 48;
		end

		if (opts.displayLevel > 1)
			fprintf('Computing gamma for problem %s.\n', name);
		end

		if (all(cellfun(@(s) isfield(s, 'gamma'), problemData)))
			if (opts.displayLevel > 1)
				fprintf(' - Already computed.\n');
			end
			continue;
		end

		[A, M, a, b, w] = loadProblem(name, problemData{1}.type);

		lambdaMin = a;
		lambdaMax = b;
		for (i = 1:1:numel(problemData))
			Di = problemData{i}.subspace.computedEigenvalues;
			lambdaMin = min(lambdaMin, min(Di));
			lambdaMax = max(lambdaMax, max(Di));
		end

		len = lambdaMax - lambdaMin;
		c = lambdaMin - 1.05*len;
		d = lambdaMax + 1.05*len;

		if (opts.parallelize)
			startParallelPool(opts.K);
		end

		[V, ~, ~, normAM, stats] = computeSubspace(A, M, c, d, opts);
		D = stats.computedEigenvalues();

		if (opts.parallelize)
			shutdownParallelPool();
		end

		for (i = 1:1:numel(problemData))
			if (isfield(problemData{i}, 'gamma'))
				if (opts.displayLevel > 1)
					fprintf(' - Already computed for run %d:  %.15e.\n', i, problemData{i}.gamma);
				end
				continue;
			end

			Di = problemData{i}.subspace.computedEigenvalues;

			indComputed = zeros(numel(Di), 1);
			indPool = 1:numel(D);
			for (j = 1:1:numel(Di))
				[~, p] = min(abs(Di(j) - D(indPool)));
				indComputed(j) = indPool(p);
				indPool(p) = [];
			end

			[maxErr, indMaxErr] = max(abs(Di - D(indComputed))./abs(D(indComputed)));
			if (maxErr > 1.0e-8)
				fprintf('WARNING:  Inaccurate eigenvalue (?) at index %d:\n', ...
					indComputed(indMaxErr));
				fprintf('val. 1 = %.15e / val. 2 = %.15e / diff = %.2e\n', ...
					Di(indMaxErr), D(indComputed(indMaxErr)), maxErr);
			end

			Duc = D;
			Duc(indComputed) = [];

			inda = Duc < a;
			indb = Duc > b;

			dist = zeros(numel(Duc), 1);
			dist(inda) = a - Duc(inda);
			dist(indb) = Duc(indb) - b;

			isInside = dist == 0;
			if (any(isInside))
				indInside = find(isInside, 1);
				lambda = Duc(indInside);
				fprintf('WARNING:  Found "uncomputed" eigenvalue in [a, b]!\n');
				fprintf('a = %.2e / b = %.2e / lambda = %.15e\n', a, b, lambda);
			end

			problemData{i}.gamma = min(dist)/(b - a);
			fprintf('j = %d:  gamma = %.2e\n', j, problemData{i}.gamma);
		end

		data(name) = {problemData};
		save(matFile, 'data');
	end
end

function printIterationExperimentsResultsTable(matFile)
	load(matFile);

	problems = keys(data);

	maxNameLength = max(cellfun(@length, problems));
	nameFormat = sprintf('%%-%ds', maxNameLength);

	fprintf([nameFormat '          Filt.     Computed                GMRES iters.         Residuals\n'], '');
	fprintf([nameFormat '     K    iters.    eigs.      gamma      Min.  Max.  Avg.     Min.      Max.\n'], 'Problem');
	fprintf('------------------------------------------------------------------------------------------------\n');
	for (n = 1:1:numel(problems))
		name = problems{n};
		problemData = data(name);
		problemData = problemData{1};

		for (i = 1:1:numel(problemData))
			if (i == round((numel(problemData) + 1)/2))
				fprintf([nameFormat '    '], name);
			else
				fprintf([nameFormat '    '], '');
			end

			ratFiltIters = problemData{i}.subspace.iter;
			numComputedEigs = length(problemData{i}.subspace.computedEigenvalues);
			minGMRESIters = min(problemData{i}.correctionIters);
			maxGMRESIters = max(problemData{i}.correctionIters);
			avgGMRESIters = mean(problemData{i}.correctionIters);
			minResidual = min(problemData{i}.residuals);
			maxResidual = max(problemData{i}.residuals);
			if (isfield(problemData{i}, 'gamma'))
				gamma = problemData{i}.gamma;
			else
				gamma = NaN;
			end
			K = problemData{i}.opts.K;
			fprintf('%2d     %2d       %4d      %.2e    %2d    %2d    %4.1f   %.1e   %.1e\n', K, ratFiltIters, numComputedEigs, gamma, minGMRESIters, maxGMRESIters, avgGMRESIters, minResidual, maxResidual);
		end

		if (n < numel(problems))
			fprintf('------------------------------------------------------------------------------------------------\n');
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function runTimingExperiments(problem, type, opts, matFile, outFile)
	maxNumCompThreads(1);

	diary(outFile);

	fprintf('---------- PROBLEM %s (%s) ----------\n', ...
		problem, char(datetime('now')));

	[A, M, wmin, wmax, w, f] = loadProblem(problem, type);
	[u, dataNewMethod] = solveShiftedSystems(A, M, f, wmin, wmax, w, opts);
	save(matFile, 'dataNewMethod');
	[u, dataNaiveMethod] = solveNaive(A, M, f, w, opts);
	save(matFile, 'dataNewMethod', 'dataNaiveMethod');

	diary('off');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function runTests()
	opts.K                       = 16;
	opts.factorization           = 'ldl';
	opts.subspace.algorithm      = 'ratfilt';
	opts.subspace.tol            = 1.0e-12;
	opts.subspace.maxIter        = 50;
	opts.subspace.insideEigsOnly = false;
	opts.solveForCorrection      = true;
	opts.solver.algorithm        = 'gmres';
	opts.solver.tol              = 1.0e-10;
	opts.solver.maxIter          = 50;
	opts.solver.restart          = [];
	opts.solver.useInitialGuess  = true;
	opts.displayLevel            = 2;
	opts.parallelize             = true;

	testStandardProblem(opts);
	testGeneralizedProblem(opts);
	testWindscreenProblem(opts);
	testBoneProblem(opts);
end

function testStandardProblem(opts)
	wmin = 0.9;
	wmax = 2.1;
	w = 1.5;

	stats = cell(4, 1);

	A = fdLaplacian1D(3^6);
	f = randn(size(A, 1), 1);
	f = f/norm(f);
	[u, stats{1}] = solveShiftedSystems(A, [], f, wmin, wmax, w, opts);

	A = fdLaplacian2D(3^3);
	f = randn(size(A, 1), 1);
	f = f/norm(f);
	[u, stats{2}] = solveShiftedSystems(A, [], f, wmin, wmax, w, opts);

	A = fdLaplacian3D(3^2);
	f = randn(size(A, 1), 1);
	f = f/norm(f);
	[u, stats{3}] = solveShiftedSystems(A, [], f, wmin, wmax, w, opts);

	w = [1.25 ; 1.5 ; 1.75];
	A = fdLaplacian1D(3^6);
	f = randn(size(A, 1), 1);
	f = f/norm(f);
	[u, stats{4}] = solveShiftedSystems(A, [], f, wmin, wmax, w, opts);
end

function testGeneralizedProblem(opts)
	wmin = 15.0;
	wmax = 200.0;
	w = 65.0;

	stats = cell(3, 1);

	[A, M] = femLaplacian2D(1.0, 1.0, 0.07);
	f = randn(size(A, 1), 1);
	f = f/norm(f);
	[u, stats{1}] = solveShiftedSystems(A, M, f, wmin, wmax, w, opts);

	[A, M] = femLaplacian3D(1.0, 1.0, 1.0, 0.25);
	f = randn(size(A, 1), 1);
	f = f/norm(f);
	[u, stats{2}] = solveShiftedSystems(A, M, f, wmin, wmax, w, opts);

	w = [40.0 ; 60.0; 80.0];
	[A, M] = femLaplacian2D(1.0, 1.0, 0.07);
	f = randn(size(A, 1), 1);
	f = f/norm(f);
	[u, stats{3}] = solveShiftedSystems(A, M, f, wmin, wmax, w, opts);
end

function testWindscreenProblem(opts)
	load('matrices/windscreen.mat');
	A = Problem.A;
	M = Problem.M;
	f = Problem.f;
	f = 0*f;
	f(1) = 1.0;
	wmin = 0.5;
	wmax = 100;
	w = linspace(wmin, wmax, 200).';
	opts.solver.tol = 1.0e-8;
	[u, stats] = solveShiftedSystems(A, M, f, wmin^2, wmax^2, w.^2, opts);
	normu = sqrt(sum(u.^2));
	semilogy(w, normu);
end

function testBoneProblem(opts)
	load('matrices/bs01.mat');
	A = Problem.A;
	M = Problem.M;
	f = Problem.f;
	wmin = 2*pi;
	wmax = 2*pi*100;
	w = linspace(wmin, wmax, 200).';
	opts.solver.tol = 1.0e-8;
	[u, stats] = solveShiftedSystems(A, M, f, wmin^2, wmax^2, w.^2, opts);
	normu = sqrt(sum(u.^2));
	semilogy(w, normu);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u, stats] = solveNaive(A, M, f, w, opts)
	standardProblem = isempty(M);
	N = size(A, 1);

	if (standardProblem)
		M = speye(size(A));
	end

	if (opts.displayLevel > 0)
		if (standardProblem)
			fprintf('Solving STANDARD shifted systems problem NAIVELY.\n');
		else
			fprintf('Solving GENERALIZED shifted systems problem NAIVELY.\n');
		end

		fprintf(' - Dimension N = %d\n', N);

		if (opts.parallelize)
			fprintf(' - Running in PARALLEL.\n');
		else
			fprintf(' - Running in SERIAL.\n');
		end
	end

	u = zeros(N, length(w));

	t = tic();

	if (opts.parallelize)
		parfor (n = 1:1:length(w))
			u(:, n) = (A - w(n)*M)\f;
		end
	else
		stats.sysTimes = zeros(length(w), 1);
		for (n = 1:1:length(w))
			tic()
			u(:, n) = (A - w(n)*M)\f;
			stats.sysTimes(n) = toc();

			if (opts.displayLevel > 0)
				fprintf('    * w = % .3e:  time:  %.2e s\n', ...
					w(n), stats.sysTimes(n));
			end
		end
	end

	stats.totalTime = toc(t);

	if (opts.displayLevel > 0)
		fprintf(' - Total time:  %2e s\n', stats.totalTime);
		fprintf(' - Finished.\n');
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u, stats] = solveShiftedSystems(A, M, f, wmin, wmax, w, opts)
	standardProblem = isempty(M);
	N = size(A, 1);

	if (opts.parallelize)
		startParallelPool(opts.K);
	end

	if (opts.displayLevel > 0)
		if (standardProblem)
			fprintf('Solving STANDARD shifted systems problem.\n');
		else
			fprintf('Solving GENERALIZED shifted systems problem.\n');
		end

		fprintf(' - Dimension N = %d\n', N);

		if (opts.parallelize)
			fprintf(' - Running in PARALLEL.\n');
		else
			fprintf(' - Running in SERIAL.\n');
		end

		switch(opts.factorization)
		case 'lu'
			fprintf(' - Using LU factorization.\n');
		case 'ldl'
			fprintf(' - Using LDL'' factorization.\n');
		otherwise
			error('Invalid factorization type.');
		end
	end

	[V, fact, xk, normAM, stats.subspace] = computeSubspace(A, M, wmin, wmax, opts);

	if (opts.parallelize)
		preconData = struct('xk', xk, 'Lk', chebfun.lagrange(xk), ...
			            'fact', fact, 'M', M, 'V', V, ...
				    'parallelize', true, ...
				    'displayLevel', opts.displayLevel);
	else
		preconData = struct('xk', xk, 'Lk', chebfun.lagrange(xk), ...
			            'fact', {fact}, 'M', M, 'V', V, ...
				    'parallelize', false, ...
				    'displayLevel', opts.displayLevel);
	end

	A1 = V'*A*V;
	if (standardProblem)
		M1 = speye(size(A1));
	else
		M1 = V'*M*V;
	end

	if (opts.displayLevel > 0)
		fprintf(' - Beginning sweep.\n');
	end

	stats.projectionTimes = zeros(length(w), 1);
	stats.correctionTimes = zeros(length(w), 1);
	stats.correctionIters = zeros(length(w), 1);
	stats.combinedTimes = zeros(length(w), 1);
	stats.residuals = zeros(length(w), 1);
	stats.totalTime = 0;

	u = zeros(N, length(w));
	for (n = 1:1:length(w))
		tic();
		f1 = V'*f;
		u1 = V*((A1 - w(n)*M1)\f1);
		t1 = toc();

		if (opts.solveForCorrection)
			if (n == 1)
				u0 = [];
			end

			t2 = tic();
			if (standardProblem)
				f2 = f - V*f1;
			else
				f2 = f - M*(V*f1);
			end
			[u2, iter] = computeCorrection(A, M, w(n), f2, V, preconData, u0, opts);
			t2 = toc(t2);

			u(:, n) = u1 + u2;

			if (opts.solver.useInitialGuess)
				u0 = u2;
			end
		else
			u(:, n) = u1;

			t2 = 0.0;
			iter = 0;
		end

		if (standardProblem)
			resid = norm(A*u(:, n) - w(n)*u(:, n) - f);
		else
			resid = norm(A*u(:, n) - w(n)*(M*u(:, n)) - f);
		end

		if (opts.displayLevel > 0)
			if (opts.solveForCorrection)
				fprintf('    * w = % .3e:  iter = %2d, resid = %.2e, time:  %.2e s + %.2e s = %.2e s\n', w(n), iter, resid, t1, t2, t1 + t2);
			else
				fprintf('    * w = % .3e:  resid = %.2e, time:  %.2e s \n', w(n), resid, t1);
			end
		end

		stats.projectionTimes(n) = t1;
		stats.correctionTimes(n) = t2;
		stats.correctionIters(n) = iter;
		stats.combinedTimes(n) = t1 + t2;
		stats.residuals(n) = resid;
		stats.totalTime = stats.totalTime + t1 + t2;
	end

	if (opts.displayLevel > 0)
		fprintf('    * Total time:  %2e s\n', stats.totalTime);
		fprintf(' - Finished.\n');
	end

	if (opts.parallelize)
		shutdownParallelPool();
	end
end

function [u, iter] = computeCorrection(A, M, w, f2, V, preconData, u0, opts)
	standardProblem = isempty(M);
	tol     = opts.solver.tol;
	maxIter = opts.solver.maxIter;
	restart = opts.solver.restart;

	fnAx = @(x) projMatvec(x, A, M, w, V);
	fnPx = @(x) precon(x, w, preconData);

	switch (opts.solver.algorithm)
	case 'minres'
		f2 = fnPx(f2);
		fnPAx = @(x) fnPx(fnAx(x));
		[u, ~, ~, iter] = minres(fnPAx, f2, tol, maxIter, [], [], u0);
	case 'gmres'
		[u, ~, ~, iter, ~] = gmres(fnAx, f2, restart, tol, maxIter, fnPx, [], u0);
		iter = iter(2);
	otherwise
		error(sprintf('Invalid value ''%s'' for opts.solver.algorithm.', opts.solver.algorithm))
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [V, fact, xk, normAM, stats] = computeSubspace(A, M, wmin, wmax, opts)
	switch (opts.subspace.algorithm)
	case 'eig'
		[V, fact, xk, normAM, stats] = computeSubspaceEig(A, M, wmin, wmax, opts);
	case 'ratfilt'
		[V, fact, xk, normAM, stats] = computeSubspaceRatfilt(A, M, wmin, wmax, opts);
	otherwise
		error(sprintf('Invalid value ''%s'' for opts.subspace.algorithm.', opts.subspace.algorithm))
	end
end

function [V, fact, xk, normAM, stats] = computeSubspaceEig(A, M, wmin, wmax, opts)
	standardProblem = isempty(M);
	K = opts.K;

	N = size(A, 1);
	if (standardProblem)
		[V, D] = eig(full(A), 'vector');
	else
		[V, D] = eig(full(A), full(M), 'vector');
	end
	ind1 = (D >= wmin) & (D <= wmax);
	V = V(:, ind1);

	normAM = max(abs(D));

	xk = chebpts(K, [wmin, wmax]);
	fact = getFactorsAtPoles(A, M, xk, opts);

	stats.iter = NaN;
	stats.computedEigenvalues = D(ind1);
end

function [V, fact, xk, normAM, stats] = computeSubspaceRatfilt(A, M, wmin, wmax, opts)
	K       = opts.K;
	tol     = opts.subspace.tol;
	maxIter = opts.subspace.maxIter;
	N = size(A, 1);

	if (isempty(M))
		M = speye(size(A));
	end

	if (opts.displayLevel > 1)
		fprintf(' - Using rational filter with K = %d poles.\n', opts.K);
	end

	if (opts.displayLevel > 1)
		fprintf(' - Estimating largest eigenvalue of (A, M).\n');
	end
	normAM = eigs(A, M, 1, 'largestabs', 'Tolerance', 1.0e-2);
	if (opts.displayLevel > 1)
		fprintf(' - Estimate:  %.16e.\n', normAM);
	end

	stats.normAM = normAM;

	if (opts.displayLevel > 1)
		fprintf(' - Counting eigenvalues (inertia).\n');
	end

	tic();
	[~, D1, ~] = ldl(A - wmax*M);
	D1 = eig(D1);
	t1 = toc();

	if (opts.displayLevel > 1)
		fprintf('    * LDL^* at w = wmax:  %.2e s\n', t1);
	end

	tic();
	[~, D2, ~] = ldl(A - wmin*M);
	D2 = eig(D2);
	t2 = toc();

	if (opts.displayLevel > 1)
		fprintf('    * LDL^* at w = wmin:  %.2e s\n', t2);
		fprintf('    * Total time:  %.2e s\n', t1 + t2);
	end

	s = full(sum(D1 < 0) - sum(D2 < 0));

	if (opts.displayLevel > 1)
		fprintf(' - Counted %d eigenvalues within interval.\n', s);
	end

	stats.eigCountTime  = t1 + t2;
	stats.eigCount = s;

	[xk, ~, vk] = chebpts(K, [wmin, wmax]);

	[fact, stats.factorization] = getFactorsAtPoles(A, M, xk, opts);

	V = randn(N, round(1.25*s));
	V = V./sqrt(sum(V.^2));

	if (opts.displayLevel > 1)
		fprintf(' - Computing subspace.\n');
	end

	stats.filterTimes = [];
	stats.svdTimes    = [];
	stats.expandTimes = [];
	stats.iterTimes   = [];
	stats.totalIterationTime = 0.0;
	stats.projectedSolveTime = 0.0;

	iter = 1;
	totalTime = 0.0;
	while (iter < maxIter)
		if (opts.displayLevel > 2)
			fprintf('    * Iter %2d:  applying filter.\n', iter);
		end

		t1 = tic();
		if (opts.parallelize)
			switch (opts.factorization)
			case 'lu'
				U = fact.U;
				L = fact.L;
				p = fact.p;
				q = fact.q;

				spmd (opts.K)
					k = labindex;
					Wk = vk(k)*lusolve(L, U, p, q, V);
					W = gplus(Wk, 1);
				end
			case 'ldl'
				L = fact.L;
				D = fact.D;
				p = fact.p;

				spmd (opts.K)
					k = labindex;
					Wk = vk(k)*ldlsolve(L, D, p, V);
					W = gplus(Wk, 1);
				end
			otherwise
				error('Invalid factorization type.');
			end

			W = W{1};
		else
			W = zeros(size(V));
			for (k = 1:1:K)
				switch (opts.factorization)
				case 'lu'
					U = fact{k}.U;
					L = fact{k}.L;
					p = fact{k}.p;
					q = fact{k}.q;

					tic();
					Wk = lusolve(L, U, p, q, V);
					t1k1 = toc();
				case 'ldl'
					L = fact{k}.L;
					D = fact{k}.D;
					p = fact{k}.p;

					tic();
					Wk = ldlsolve(L, D, p, V);
					t1k1 = toc();
				otherwise
					error('Invalid factorization type.');
				end

				tic();
				W = W + vk(k)*Wk;
				t1k2 = toc();

				if (opts.displayLevel > 2)
					fprintf('       + Pole %2d:  solve:  %.2e s / accumulate:  %.2e s\n', k, t1k1, t1k2);
				end
			end
		end
		t1 = toc(t1);

		tic();
		[V, S, ~] = svd(W, 0);
		S = diag(S);
		t2 = toc();

		totalTime = totalTime + t1 + t2;

		if (opts.displayLevel > 1)
			fprintf('    * Iter %2d:  dim = %d, ratio = %.2e / filter:  %.2e s / SVD:  %.2e s\n', iter, size(V, 2), S(end)/S(1), t1, t2);
		end

		if (S(end) < tol*S(1))
			tic();

			A1 = V'*A*V;
			A1 = (A1 + A1')/2.0;
			M1 = V'*M*V;
			M1 = (M1 + M1')/2.0;
			[V1, D1] = eig(A1, M1, 'vector');
			VV = V*V1;
			D = D1;

			R = A*VV - M*(VV.*(D.'));
			resid = sqrt(sum(R.^2));
			ind = resid < tol*normAM.*sqrt(sum(VV.^2));
			VV = VV(:, ind);
			D = D(ind);

			t = toc();

			if (opts.displayLevel > 1)
				fprintf('    * Compute eigenvectors / residual check:  %.2e s.\n', t);
			end

			stats.projectedSolveTime = stats.projectedSolveTime + t;

			numEigsInInterval = sum((D >= wmin) & (D <= wmax));
			if (numEigsInInterval ~= s)
				fprintf('    * Found %d eigenvalues in interval; expected %d.  Continuing to iterate.\n', numEigsInInterval, s);
			else
				stats.filterTimes = [stats.filterTimes ; t1];
				stats.svdTimes    = [stats.svdTimes ; t2];
				stats.expandTimes = [stats.expandTimes;  0.0];
				stats.iterTimes   = [stats.iterTimes ; t1 + t2 + 0.0];
				stats.totalIterationTime = totalTime;
				stats.iter = iter;
				stats.totalTime = stats.totalIterationTime + stats.projectedSolveTime;

				if (opts.displayLevel > 1)
					fprintf('    * Total time:  %.2e s\n', stats.totalTime);
					fprintf(' - After checking residuals, %d eigenvectors remain.\n', size(VV, 2));
				end

				V = VV;

				if (opts.subspace.insideEigsOnly)
					if (opts.displayLevel > 1)
						fprintf(' - Keeping %d eigenvectors from inside the interval only.\n', s);
				end
					ind = (D >= wmin) & (D <= wmax);
					D = D(ind);
					V = V(:, ind);
				end

				stats.computedEigenvalues = D;

				break;
			end
		end

		tic();
		v = randn(N, round(0.1*size(V, 2)));
		V = [V v/norm(v)];
		t3 = toc();

		totalTime = totalTime + t3;

		stats.filterTimes = [stats.filterTimes ; t1];
		stats.svdTimes    = [stats.svdTimes ; t2];
		stats.expandTimes = [stats.expandTimes;  t3];
		stats.iterTimes   = [stats.iterTimes ; t1 + t2 + t3];

		iter = iter + 1;
	end
end

function [fact, stats] = getFactorsAtPoles(A, M, xk, opts)
	if (opts.parallelize)
		[fact, stats] = getFactorsAtPolesParallel(A, M, xk, opts);
	else
		[fact, stats] = getFactorsAtPolesSerial(A, M, xk, opts);
	end
end

function [fact, stats] = getFactorsAtPolesSerial(A, M, xk, opts)
	K = length(xk);

	if (opts.displayLevel > 1)
		fprintf(' - Factoring pencil at filter poles (serial).\n');
	end

	stats.factTime = zeros(K, 1);
	stats.factMem = zeros(K, 1);
	stats.totalTime = 0;
	stats.totalMemory = 0;

	fact = cell(K, 1);
	for (k = 1:1:K)
		tic();

		switch (opts.factorization)
		case 'lu'
			[L, U, p, q] = lu(A - xk(k)*M, 'vector');
			fact{k} = struct('L', L, 'U', U, 'p', p, 'q', q, 'type', 'lu');
		case 'ldl'
			[L, D, p] = ldl(A - xk(k)*M, 'vector');
			fact{k} = struct('L', L, 'D', D, 'p', p, 'type', 'ldl');
		otherwise
			error('Invalid factorization type.');
		end

		t = toc();

		stats.totalTime = stats.totalTime + t;
		if (opts.displayLevel > 1)
			switch (opts.factorization)
			case 'lu'
				memL = whos('L').bytes;
				memU = whos('U').bytes;
				memp = whos('p').bytes;
				memq = whos('q').bytes;
				mem = memL + memU + memp + memq;
			case 'ldl'
				memL = whos('L').bytes;
				memD = whos('D').bytes;
				memp = whos('p').bytes;
				mem = memL + memD + memp;
			otherwise
				error('Invalid factorization type.');
			end

			digits = floor(log10(K)) + 1;
			fprintf(sprintf('    * Fact. at pole %%%dd:  %%.2e s / %%s\n', digits), k, t, formatMemory('%.2f', mem));

			stats.factTime(k) = t;
			stats.factMem(k) = mem;
			stats.totalTime = stats.totalTime + t;
			stats.totalMemory = stats.totalMemory + mem;
		end
	end

	if (opts.displayLevel > 1)
		fprintf('    * Total time/memory:  %.2e s / %s\n', stats.totalTime, formatMemory('%.2f', stats.totalMemory));
	end

	stats.totalTime = stats.totalTime;
	stats.totalMemory = stats.totalMemory;
end

function [fact, stats] = getFactorsAtPolesParallel(A, M, xk, opts)
	K = length(xk);

	if (opts.displayLevel > 1)
		fprintf(' - Factoring pencil at filter poles (parallel).\n');
	end

	switch (opts.factorization)
	case 'lu'
		tic();

		spmd (opts.K)
			k = labindex;
			[L, U, p, q] = lu(A - xk(k)*M, 'vector');
		end

		totalTime = toc();

		if (opts.displayLevel > 1)
			fprintf('    * Total time:  %.2e s\n', totalTime);
		end

		fact = struct('L', L, 'U', U, 'p', p, 'q', q, 'type', 'lu');
	case 'ldl'
		tic();

		spmd (opts.K)
			k = labindex;
			[L, D, p] = ldl(A - xk(k)*M, 'vector');
		end

		totalTime = toc();

		if (opts.displayLevel > 1)
			fprintf('    * Total time:  %.2e s\n', totalTime);
		end

		fact = struct('L', L, 'D', D, 'p', p, 'type', 'ldl');
	otherwise
		error('Invalid factorization type.');
	end

	stats.totalTime = totalTime;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function y = projMatvec(x, A, M, w, V)
	if (isempty(M))
		y = x - V*(V'*x);
		y = A*y - w*y;
		y = y - V*(V'*y);
	else
		y = x - V*(V'*(M*x));
		y = A*y - w*(M*y);
		y = y - M*(V*(V'*y));
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Y = precon(X, w, preconData)
	xk   = preconData.xk;
	fact = preconData.fact;
	M    = preconData.M;
	V    = preconData.V;
	K    = length(xk);

	standardProblem = isempty(M);

	if (preconData.displayLevel > 2)
		fprintf('    * Applying preconditioner.\n');
	end

	ck = feval(preconData.Lk, w);

	t = tic();
	if (standardProblem)
		X = X - V*(V'*X);
	else
		X = X - M*(V*(V'*X));
	end
	t = toc(t);

	if (preconData.displayLevel > 2)
		fprintf('       + Initial projection:  %.2e s\n', t);
	end

	if (preconData.parallelize)
		switch (fact.type)
		case 'lu'
			U = fact.U;
			L = fact.L;
			p = fact.p;
			q = fact.q;

			spmd (K)
				k = labindex;
				Yk = ck(k)*lusolve(L, U, p, q, X);
				Y = gplus(Yk, 1);
			end
			Y = Y{1};
		case 'ldl'
			L = fact.L;
			D = fact.D;
			p = fact.p;

			spmd (K)
				k = labindex;
				Yk = ck(k)*ldlsolve(L, D, p, X);
				Y = gplus(Yk, 1);
			end
			Y = Y{1};
		otherwise
			error('Invalid factorization type.');
		end
	else
		Y = zeros(size(X));
		for (k = 1:1:K)
			switch (fact{k}.type)
			case 'lu'
				U = fact{k}.U;
				L = fact{k}.L;
				p = fact{k}.p;
				q = fact{k}.q;

				t1 = tic();
				Yk = lusolve(L, U, p, q, X);
				t1 = toc(t1);
			case 'ldl'
				L = fact{k}.L;
				D = fact{k}.D;
				p = fact{k}.p;

				t1 = tic();
				Yk = ldlsolve(L, D, p, X);
				t1 = toc(t1);
			otherwise
				error('Invalid factorization type.');
			end

			t2 = tic();
			Y = Y + ck(k)*Yk;
			t2 = toc(t2);

			if (preconData.displayLevel > 2)
				fprintf('       + Pole %2d:  solve:  %.2e s / accumulate:  %.2e s\n', k, t1, t2);
			end
		end
	end

	t = tic();
	if (standardProblem)
		Y = Y - V*(V'*Y);
	else
		Y = Y - V*(V'*(M*Y));
	end
	t = toc(t);

	if (preconData.displayLevel > 2)
		fprintf('       + Terminal projection:  %.2e s\n', t);
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X = lusolve(L, U, p, q, B)
	X = U\(L\B(p, :));
	X(q, :) = X;
end

function X = ldlsolve(L, D, p, B)
	X = L\B(p, :);
	d = spparms('bandden');
	spparms('bandden', 0);
	X = D\X;
	spparms('bandden', d);
	X = L'\X;
	X(p, :) = X;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function L = fdLaplacian1D(n)
	e = ones(n, 1);
	L = spdiags([-e 2*e -e], [-1 0 1], n, n);
end

function L = fdLaplacian2D(n)
        I = speye(n);
        L1 = fdLaplacian1D(n);
        L = kron(L1, I) + kron(I, L1);
end

function L = fdLaplacian3D(n)
        L1 = fdLaplacian1D(n);
        L2 = fdLaplacian2D(n);
        L = kron(L2, speye(n)) + kron(speye(n^2), L1);
end

function [A, M] = femLaplacian2D(L, W, Hmax)
	vx = [0 ; L ; L ; 0];
	vy = [0 ; 0 ; W ; W];
	T = [2 ; 4 ; vx ; vy];
	geom = decsg(T);

	model = createpde();
	geometryFromEdges(model, geom);
	generateMesh(model, 'Hmax', Hmax);
	specifyCoefficients(model, 'm', 0, 'c', 1, 'a', 0, 'd', 1, 'f', 0);

	FEM = assembleFEMatrices(model);
	A = FEM.K;
	M = FEM.M;

	A = (A + A')/2;
	M = (M + M')/2;
end

function [A, M] = femLaplacian3D(L, W, H, Hmax)
	model = createpde();
	model.Geometry = multicuboid(L, W, H);
	generateMesh(model, 'Hmax', Hmax);
	specifyCoefficients(model, 'm', 0, 'c', 1, 'a', 0, 'd', 1, 'f', 0);

	FEM = assembleFEMatrices(model);
	A = FEM.K;
	M = FEM.M;

	A = (A + A')/2;
	M = (M + M')/2;
end

function findFEMMeshParameters()
	n = [1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6];

	hMin2D = 0.002;
	hMax2D = 0.03;
	femFunction2D = @(h) femLaplacian2D(1.0, 2.0^0.25, h);
	for (i = 1:1:length(n))
		fprintf(' - Finding h for 2D, n = %d.\n', n(i));
		findFEMMeshParameter(n(i), hMin2D, hMax2D, femFunction2D);
	end

	hMin3D = 0.02;
	hMax3D = 0.2;
	femFunction3D = @(h) femLaplacian3D(1.0, 2.0^0.25, 3.0^0.25, h);
	for (i = 1:1:length(n))
		fprintf(' - Finding h for 3D, n = %d.\n', n(i));
		findFEMMeshParameter(n(i), hMin3D, hMax3D, femFunction3D);
	end
end

function h = findFEMMeshParameter(nTarget, hMin, hMax, femLaplacian)
	if (nTarget < 1e5)
		tol = 1.0e-2;
	else
		tol = 1.0e-3;
	end

	while (1)
		h = (hMin + hMax)/2.0;
		fprintf('    * Trying h = %.15e... ', h);
		[A, M] = femLaplacian(h);
		n = size(A, 1);
		fprintf('dimension %d\n', n);

		if ((abs(n - nTarget)/nTarget < tol) || (abs(n - nTarget) < 5.0*10*tol))
			fprintf('    * Success with h = %.15e.\n', h);
			break;
		elseif (n > nTarget)
			hMin = h;
		elseif (n < nTarget)
			hMax = h;
		end
	end
end

function generateFEMMatrices(regenerate)
	if (nargin < 1)
		regenerate = false;
	end

	femLaplacian2DFunc = @(Hmax) femLaplacian2D(1.0, 2.0^0.25, Hmax);
	femLaplacian3DFunc = @(Hmax) femLaplacian3D(1.0, 2.0^0.25, 3.0^0.25, Hmax);

	femProblems = {...
	        'fem2D-10K',  struct('func', femLaplacian2DFunc, 'Hmax', 2.354687500000000e-02), ...
	        'fem2D-20K',  struct('func', femLaplacian2DFunc, 'Hmax', 1.665625000000000e-02), ...
		'fem2D-50K',  struct('func', femLaplacian2DFunc, 'Hmax', 1.053125000000000e-02), ...
		'fem2D-100K', struct('func', femLaplacian2DFunc, 'Hmax', 7.409426078490841e-03), ...
		'fem2D-200K', struct('func', femLaplacian2DFunc, 'Hmax', 5.199218750000002e-03), ...
		'fem2D-500K', struct('func', femLaplacian2DFunc, 'Hmax', 3.252685546875000e-03), ...
		'fem2D-1M',   struct('func', femLaplacian2DFunc, 'Hmax', 2.283691406250001e-03), ...
		'fem3D-10K',  struct('func', femLaplacian3DFunc, 'Hmax', 1.240625000000000e-01), ...
		'fem3D-20K',  struct('func', femLaplacian3DFunc, 'Hmax', 9.664062500000001e-02), ...
		'fem3D-50K',  struct('func', femLaplacian3DFunc, 'Hmax', 7.132812499999999e-02), ...
		'fem3D-100K', struct('func', femLaplacian3DFunc, 'Hmax', 5.625488281250000e-02), ...
		'fem3D-200K', struct('func', femLaplacian3DFunc, 'Hmax', 4.461908781234338e-02), ...
		'fem3D-500K', struct('func', femLaplacian3DFunc, 'Hmax', 3.279128313064574e-02), ...
		'fem3D-1M',   struct('func', femLaplacian3DFunc, 'Hmax', 2.605895996093750e-02), ...
	};
	femProblems = dictionary(femProblems{:});

	for (name = keys(femProblems).')
		fprintf('Generating matrix %s... ', name);
		matFile = sprintf('matrices/%s.mat', name);
		if ((exist(matFile) == 2) && ~regenerate)
			fprintf('already exists.\n');
			continue;
		end

		femLaplacianFunc = femProblems(name).func;
		Hmax = femProblems(name).Hmax;

		tic();
		[A, M] = femLaplacianFunc(Hmax);
		t = toc();

		Problem.A = A;
		Problem.M = M;
		save(matFile, 'Problem');

		fprintf('Done.  (time = %.2e s, n = %d)\n', t, size(A, 1));
	end
end

function [A, M, wmin, wmax, w, f] = loadProblem(name, type)
	switch (name)
	case 'Kuu/Muu'
		load('matrices/Kuu.mat');
		A = Problem.A;
		load('matrices/Muu.mat');
		M = Problem.A;
		switch (type)
		case 'low'
			wmin = 9.294244596423432e+00;
			wmax = 1.671034707157513e+03;
		case 'mid'
			wmin = 7.5e3;
			wmax = 9.5e3;
		end
		w = linspace(wmin, wmax, 100).';
	case 'crystk03/crystm03'
		load('matrices/crystk03.mat');
		A = Problem.A;
		load('matrices/crystm03.mat');
		M = Problem.A;
		c = norm(A, 'fro');
		A = A/c;
		M = M/c;
		switch (type)
		case 'low'
			wmin = -1.0e-01;
			wmax =  7.152557373046977e-01;
		case 'mid'
			wmin = 1.0e3;
			wmax = 1.25e3;
		end
		w = linspace(wmin, wmax, 100).';
	case 'bcsstk39/bcsstm39'
		load('matrices/bcsstk39.mat');
		A = Problem.A;
		load('matrices/bcsstm39.mat');
		M = Problem.A;
		switch (type)
		case 'low'
			error('bcsstk39/bcsstm39 is not positive definite.');
		case 'mid'
			wmin = 1.0e5;
			wmax = 1.25e5;
		end
		w = linspace(wmin, wmax, 100).';
	case 'qa8fk/qa8fm'
		load('matrices/qa8fk.mat');
		A = Problem.A;
		load('matrices/qa8fm.mat');
		M = Problem.A;
		switch (type)
		case 'low'
			wmin = -1.0e-1;
			wmax =  3.285705566406255e+01;
		case 'mid'
			wmin = 1.0e3;
			wmax = 1.01e3;
		end
		w = linspace(wmin, wmax, 100).';
	case 'gyro_k/gyro_m'
		load('matrices/gyro_k.mat');
		A = Problem.A;
		load('matrices/gyro_m.mat');
		M = Problem.A;
		switch (type)
		case 'low'
			error('Type low not supported for gyro_k/gyro_m.');
		case 'mid'
			wmin = 1.0e15;
			wmax = 2.0e15;
		end
		w = linspace(wmin, wmax, 100).';
	case 'windscreen'
		load('matrices/windscreen.mat');
		A = Problem.A;
		M = Problem.M;
		switch (type)
		case 'low'
			wmin = -1.0e-1;
			wmax =  9.130859375000149e+06;
		case 'mid'
			wmin = 1.0e7;
			wmax = 5.0e7;
		end
		w = linspace(wmin, wmax, 100).';
	case 'bs01'
		load('matrices/bs01.mat');
		A = Problem.A;
		M = Problem.M;
		switch (type)
		case 'low'
			wmin = 2.963469920052067e+03;
			wmax = 5.226602711557303e+07;
		case 'mid'
			wmin = 1.0e7;
			wmax = 6.0e7;
		end
		w = linspace(wmin, wmax, 100).';
	case 'b010'
		load('matrices/b010.mat');
		A = Problem.A;
		M = Problem.M;
	case 'fem2D-10K'
		load('matrices/fem2D-10K.mat');
		A = Problem.A;
		M = Problem.M;
		switch (type)
		case 'low'
			wmin = -1.0e-1;
			wmax =  1.801757812500000e+03;
		case 'mid'
			error('Mid-frequency for FEM matrices not yet supported.');
		end
		w = linspace(wmin, wmax, 100).';
	case 'fem2D-10K-T'
		load('matrices/fem2D-10K.mat');
		A = Problem.A;
		M = Problem.M;
		switch (type)
		case 'low'
			wmin = -1.0e-1;
			wmax =  1.801757812500000e+03;
		case 'mid'
			error('Mid-frequency for FEM matrices not yet supported.');
		end
		w = linspace(wmin, wmax, 10).';
	case 'fem2D-20K'
		load('matrices/fem2D-20K.mat');
		A = Problem.A;
		M = Problem.M;
		switch (type)
		case 'low'
			wmin = -1.0e-1;
			wmax =  1.796875000000000e+03;
		case 'mid'
			error('Mid-frequency for FEM matrices not yet supported.');
		end
		w = linspace(wmin, wmax, 100).';
	case 'fem2D-50K'
		load('matrices/fem2D-50K.mat');
		A = Problem.A;
		M = Problem.M;
		switch (type)
		case 'low'
			wmin = -1.0e-1;
			wmax =  1.796875000000000e+03;
		case 'mid'
			error('Mid-frequency for FEM matrices not yet supported.');
		end
		w = linspace(wmin, wmax, 100).';
	case 'fem2D-100K'
		load('matrices/fem2D-100K.mat');
		A = Problem.A;
		M = Problem.M;
		switch (type)
		case 'low'
			wmin = -1.0e-1;
			wmax =  1.796875000000000e+03;
		case 'mid'
			error('Mid-frequency for FEM matrices not yet supported.');
		end
		w = linspace(wmin, wmax, 100).';
	case 'fem2D-200K'
		load('matrices/fem2D-200K.mat');
		A = Problem.A;
		M = Problem.M;
		switch (type)
		case 'low'
			wmin = -1.0e-1;
			wmax =  1.796875000000000e+03;
		case 'mid'
			error('Mid-frequency for FEM matrices not yet supported.');
		end
		w = linspace(wmin, wmax, 100).';
	case 'fem2D-500K'
		load('matrices/fem2D-500K.mat');
		A = Problem.A;
		M = Problem.M;
		switch (type)
		case 'low'
			wmin = -1.0e-1;
			wmax =  1.796875000000000e+03;
		case 'mid'
			error('Mid-frequency for FEM matrices not yet supported.');
		end
		w = linspace(wmin, wmax, 100).';
	case 'fem2D-1M'
		load('matrices/fem2D-1M.mat');
		A = Problem.A;
		M = Problem.M;
		switch (type)
		case 'low'
			wmin = -1.0e-1;
			wmax =  1.796875000000000e+03;
		case 'mid'
			error('Mid-frequency for FEM matrices not yet supported.');
		end
		w = linspace(wmin, wmax, 100).';
	case 'fem3D-10K'
		load('matrices/fem3D-10K.mat');
		A = Problem.A;
		M = Problem.M;
		switch (type)
		case 'low'
			wmin = -1.0e-1;
			wmax =  2.968750000000000e+02;
		case 'mid'
			error('Mid-frequency for FEM matrices not yet supported.');
		end
		w = linspace(wmin, wmax, 100).';
	case 'fem3D-20K'
		load('matrices/fem3D-20K.mat');
		A = Problem.A;
		M = Problem.M;
		switch (type)
		case 'low'
			wmin = -1.0e-1;
			wmax =  2.949218750000000e+02;
		case 'mid'
			error('Mid-frequency for FEM matrices not yet supported.');
		end
		w = linspace(wmin, wmax, 100).';
	case 'fem3D-50K'
		load('matrices/fem3D-50K.mat');
		A = Problem.A;
		M = Problem.M;
		switch (type)
		case 'low'
			wmin = -1.0e-1;
			wmax =  2.929687500000000e+02;
		case 'mid'
			error('Mid-frequency for FEM matrices not yet supported.');
		end
		w = linspace(wmin, wmax, 100).';
	case 'fem3D-100K'
		load('matrices/fem3D-100K.mat');
		A = Problem.A;
		M = Problem.M;
		switch (type)
		case 'low'
			wmin = -1.0e-1;
			wmax =  2.929687500000000e+02;
		case 'mid'
			error('Mid-frequency for FEM matrices not yet supported.');
		end
		w = linspace(wmin, wmax, 100).';
	case 'fem3D-200K'
		load('matrices/fem3D-200K.mat');
		A = Problem.A;
		M = Problem.M;
		switch (type)
		case 'low'
			wmin = -1.0e-1;
			wmax =  2.929687500000000e+02;
		case 'mid'
			error('Mid-frequency for FEM matrices not yet supported.');
		end
		w = linspace(wmin, wmax, 100).';
	case 'fem3D-500K'
		load('matrices/fem3D-500K.mat');
		A = Problem.A;
		M = Problem.M;
		switch (type)
		case 'low'
			wmin = -1.0e-1;
			wmax =  2.929687500000000e+02;
		case 'mid'
			error('Mid-frequency for FEM matrices not yet supported.');
		end
		w = linspace(wmin, wmax, 100).';
	case 'fem3D-500K-T'
		load('matrices/fem3D-500K.mat');
		A = Problem.A;
		M = Problem.M;
		switch (type)
		case 'low'
			wmin = -1.0e-1;
			wmax =  2.929687500000000e+02;
		case 'mid'
			error('Mid-frequency for FEM matrices not yet supported.');
		end
		w = linspace(wmin, wmax, 1000).';
	case 'fem3D-1M'
		load('matrices/fem3D-1M.mat');
		A = Problem.A;
		M = Problem.M;
		switch (type)
		case 'low'
			wmin = -1.0e-1;
			wmax =  2.929687500000000e+02;
		case 'mid'
			error('Mid-frequency for FEM matrices not yet supported.');
		end
		w = linspace(wmin, wmax, 100).';
	otherwise
		error(sprintf('Unknown problem ''%s''.\n', name))
	end

	f = randn(size(A, 1), 1);
	f = f/norm(f);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function c = findEigenvalueCutoff(A, M, n, cMax)
	cMin = 0.0;

	if (nargin < 4)
		fprintf('    * Computing crude upper bound... ');

		tic();
		%cMax = eigs(A, M, 1, 'largestabs', 'Tolerance', 1.0e-2);
		cMax = normest(A, 1e-3)*condest(M)/normest(M, 1e-3);
		t = toc();

		fprintf('cMax = %.15e (%.2e s)\n', cMax, t);
	end

	while (1)
		c = (cMin + cMax)/2.0;
		fprintf('    * Trying c = %.15e... ', c);

		tic()
		[~, D1, ~] = ldl(A - c*M);
		D1 = eig(D1);
		t = toc();

		nc = full(sum(D1 < 0));
		fprintf('%d eigenvalues (%.2e s)\n', nc, t);

		if (nc > n)
			cMax = c;
		elseif (nc < n)
			cMin = c;
		else
			fprintf('    * Found c = %.15e\n', c);
			break;
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function startParallelPool(K)
	pool = gcp('nocreate');
	if (isempty(pool) || (pool.NumWorkers ~= K))
		delete(pool);
		pc = parcluster('local');
		pc.NumWorkers = K;
		parpool('local', K, 'IdleTimeout', 360);
	end
end

function shutdownParallelPool(K)
	pool = gcp('nocreate');
	delete(pool);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s = formatMemory(floatFormat, b)
	prefixes = {'', 'K', 'M', 'G', 'T'};

	n = 1;
	while (b > 1024)
		b = b/1024;
		n = n + 1;
	end

	s = sprintf(floatFormat, b);
	s = sprintf('%s %sB', s, prefixes{n});
end
