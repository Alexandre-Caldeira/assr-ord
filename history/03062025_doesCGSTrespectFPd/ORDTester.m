classdef ORDTester 
    properties
        % Test Parameters     
        NDC 
        desired_alpha
        

        % Exam/ORD data
        epochs = []
        epochs_index_metadata
        epochs_method
        K_stages    % total number of tests to be applied on the ORD
        nWindows    % number of Windows (may match K, or not)

        % Esses já existem no ORDCalculator: 
        % (deveria estar aqui?)
        subjectIndices = [1:11];
        stimulusIndices = [1:5];
        startWindows = [1];
        windowStepSizes = [24 32];
        lastWindows = [50];

        % (bad practice, fix later)
        dataloader
        ord_calculator

        % Utils
        timer   % obj timer
        id      % obj id
    end


    methods
        function obj = ORDTester(ord_calculator,p)
            arguments
                ord_calculator

                % Optional inputs
                p.dataloader = [];
                p.K_stages = [];
                p.desired_alpha = 0.05;              

            end

            obj.timer = tic;
            obj.id = [num2str(keyHash(obj.timer))];
            obj.ord_calculator = ord_calculator;

            if isempty(p.K_stages)
                p.K_stages = ord_calculator.K_stages;
            end

            if isa(p.dataloader, 'double')
                obj.dataloader = obj.ord_calculator.dataloader;
            else
                obj.dataloader = p.dataloader;
            end

        end

        function obj = compute_beta_cgst_thresholds(obj)
            warning('Test pending')

            

        end

        function obj = compute_chesnaye_thresholds(obj)
            error('Implementation pending')
        end

        function obj = compute_zanoteli_thresholds(obj)
            error('Implementation pending')
        end

        % UTILS
        function age(obj)
            fprintf( ...
                '\n\tThis ORDTester was built %0.2f seconds ago.\n\n', ...
                round(toc(obj.timer),2))
        end


    end
end