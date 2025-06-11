classdef disper

    properties
        id;
        text;
        timer;
    end
    
    methods

        function obj = disper(text)
            obj.timer = tic;
            obj.text = text;
            obj.id = [num2str(keyHash(keyHash(text)+obj.timer))];
        end
        
        function show(obj)
            fprintf('\n%s\n\n',obj.text)
        end
        
        function age(obj)
            fprintf('\n%f secs\n\n', toc(obj.timer))
        end

        function disp(obj, varargin)

            if nargin==1
                age_now = toc(obj.timer);
                fprintf( ...
                    'This disper is %f secs old, and says: \n%s\n\n', ...
                    age_now,obj.text);

            elseif nargin==2
                friend = varargin{1};
                if ~isa(friend, 'disper')
                    error('Friend must be a disper object');
                end

                fprintf('This disper is: \n\t> %s\n\n and %s is a friend!\n\n', ...
                    obj.id, friend.id)

                disp(obj)

                disp(friend)

            else 
                fprintf('That is just too many friends!\n')
            end
        end

    end

end