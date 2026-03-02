classdef Node < handle

    properties
        rank
        children
        no_of_children
        eigenvalues
        V_rows
        U_columns
        rowsWithZero
    end
    
    methods
        function obj = Node()
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.no_of_children = 0;
        end
        
        function obj = append_child(obj,node)
            if obj.no_of_children == 0
                obj.children = [node];
                obj.no_of_children = obj.no_of_children + 1;
            else
                obj.children(end+1) = node;
                obj.no_of_children = obj.no_of_children + 1;
            end

        end
    end
end



