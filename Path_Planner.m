clc;
clear all;
%% Plotting Grid
hold on;
Grid = ones(11,11);
Grid(2,1) = 0;
Grid(1:2,4) = 0;
Grid(4,3:4) = 0;
Grid(5,4) = 0;
Grid(7,5) = 0;
Grid(6,9) = 0; %Grid(6,7) = 0;
Grid(8:9,7:8) = 0;
%Grid(9:10,2)=0;
%Grid(10,6) = 0;
Grid(2:3,7) = 0; %Grid(2:3,7) = 0;
%Grid(1:3,10)=0;
%Grid(1,9) = 0;
%Grid(5,9) = 0;
pcolor(Grid)
colormap(gray(2))
% To flip the vertical axes, uncomment next line
% axis ij;
text(1.25,1.5,'Start')
text(10.25,10.5,'Goal')

% Input locations
Set_X = 2.5:1:9.5;
Set_Y = 1.5:1:9.5;
Chromosome_Length = 4;% Excluding first and last points
Population_Size = 20;% Mating population size
Initial_Population = cell(Population_Size,1);
Crossover_Population = cell(Population_Size,1);
Offspring_Population = cell(Population_Size,1);
Feasible_Offspring = cell(Population_Size,1);
epoch = 12;

Best_Generation = cell(epoch,1);
fitness_best = [];
best_population = 0;
crossover_probability = 1;
mutation_probability = 0.9;
Total_Population = cell(2*Population_Size,1);
% Initial Population Generation
i=1;c=1;collect_distinguish = {};
while i<Population_Size+1
    Population = [];
    id1 = randi([1 8],1,Chromosome_Length);
    id2 = randi([1 9],1,Chromosome_Length);
    X = Set_X(id1);
    Y = Set_Y(id2);
    Chromosome =[];
    for j = 1:1:Chromosome_Length
        Chromosome = [Chromosome;X(j) Y(j)];
    end
    % Population = [1.5 2.5 3.5 3.5 5.5 7.5 10.5 10.5; 1.5 1.5 2.5 3.5 3.5 5.5 5.5 10.5];
    Population = ([[1.5 1.5];Chromosome;[10.5 10.5]])';
    inarg = [Population(1,:);Population(2,:)];
    Distinguishing_Factor = distinguish_algo_v5(Population);
    if sum(Distinguishing_Factor(1,Chromosome_Length+3:end))==0
        Initial_Population{i}=fitness(Population);
        i = i+1;
    end
    collect_distinguish{c} = Distinguishing_Factor;
    c=c+1;
    
end

%end
%% GENETIC ALGORITHM STARTS
for generation = 1:1:epoch
    index = 1;
    fitness_values = [];
   if generation==1 || generation==2 
       max_num = 0;
   else 
        max_num = 6;
   end
   
    for i = 1:1:Population_Size
        Crossover_Population{i,1} = Initial_Population{i,1}(:,1:end-1);
    end
    while index<Population_Size+1
        for i = 1:2:Population_Size-1
            A = Crossover_Population{i,1};
            B = Crossover_Population{i+1,1};
            if rand(1)< crossover_probability
                random_id = randi([200 100*(Chromosome_Length+1)],1,1);
                random_id = fix((random_id)/100);
                temp = A(:,random_id);
                A(:,random_id) = B(:,random_id);
                B(:,random_id) = temp;
            end
            if rand(1)< mutation_probability
                random_id1 = randi([200 100*(Chromosome_Length+1)],1,1);
                random_id1 = fix((random_id1)/100);
                random_id2 = randi([200 100*(Chromosome_Length+1)],1,1);
                random_id2 = fix((random_id2)/100);
                temp = A(:,random_id1);
                A(:,random_id1) = A(:,random_id2);
                A(:,random_id2) = temp;
            end
            if rand(1)< mutation_probability
                random_id1 = randi([2 Chromosome_Length+1],1,1);
                random_id2 = randi([2 Chromosome_Length+1],1,1);
                temp = B(:,random_id1);
                B(:,random_id1) = B(:,random_id2);
                B(:,random_id2) = temp;
            end
            Offspring_Population{i,1} = A;
            Offspring_Population{i+1,1} = B;
        end
        for i=1:1:Population_Size
            Distinguishing_Factor = distinguish_algo_v5(Offspring_Population{i,1});
            if sum(Distinguishing_Factor(1,Chromosome_Length+3:end))==0
                Feasible_Offspring{index} = fitness(Offspring_Population{i,1});
                index = index+1;
            end
            
        end
    end
    
    Total_Population = {[Initial_Population;Feasible_Offspring]};
    Total_Population = Total_Population{1,1};
    for i = 1:1:length(Total_Population)
        fitness_values(end+1) = Total_Population{i,1}(1,end);
    end
    [fitness1 id1_fitness] = maxk(fitness_values,max_num);
    [fitness2 id2_fitness] = mink(fitness_values,Population_Size-max_num);
    id_fitness = [id1_fitness id2_fitness];
    for i=1:1:Population_Size
        Initial_Population{i,1} = Total_Population{id_fitness(i),1}(:,1:end);
    end
    Best_Generation{generation} = Initial_Population{1,1};
    fitness_best(generation) = Initial_Population{1,1}(1,end);
    
    % Schuffling the population now
    change_order = randperm(Population_Size);
    New_Population_Schuffled = cell(Population_Size,1); local = 0;
    for j = change_order
        local = local+1;
        New_Population_Schuffled{local,1} = Initial_Population{j,1};
    end
    Initial_Population = New_Population_Schuffled;
    
end
hold on
plot(Best_Generation{10,1}(1,1:end-1),Best_Generation{10,1}(2,1:end-1));
figure(2)
hold off
plot(1:1:epoch,fitness_best);
title('fitness vs generation');
xlabel('Generation');
ylabel('Fitness');

%% DA function
function f = distinguish_algo_v5(chromosome) % path as argument

%% Obstacle border line equations
O1 = [0 1 2; 0 1 3; 1 0 1; 1 0 2]; % col 1 row 2
O2 = [0 1 1; 0 1 3; 1 0 4; 1 0 5];  % col 4 row 1:2
O3_1 = [0 1 4; 0 1 5; 1 0 3; 1 0 5]; % col 3:4 row 4
O3_2 = [0 1 5; 0 1 6; 1 0 4; 1 0 5]; %col 4 row 5
O4 = [0 1 2; 0 1 4; 1 0 7; 1 0 8]; %col 7 row 2:3
O5 = [0 1 6; 0 1 7; 1 0 9; 1 0 10]; %col 7 row 6  [0 1 6; 0 1 7; 1 0 7; 1 0 8];
O6 = [0 1 7; 0 1 8; 1 0 5; 1 0 6]; %col 5 row 7
O7 = [0 1 8; 0 1 10; 1 0 7; 1 0 9]; % col 7:8 row 8:9
obstacle_xy = [O1(1:4,1:2); O2(1:4,1:2); O3_1(1:4,1:2); O3_2(1:4,1:2); O4(1:4,1:2);...
    O5(1:4,1:2);  O6(1:4,1:2);  O7(1:4,1:2)]';
obstacle_c =  [O1(1:4,3); O2(1:4,3); O3_1(1:4,3); O3_2(1:4,3); O4(1:4,3);...
    O5(1:4,3);  O6(1:4,3);  O7(1:4,3)]';
for i = 1:1:length(chromosome)-1    %generating line segment using formula
    count = 1;
    for j = 1:1:length(obstacle_c)
        xy_1 = chromosome(1:2,i);       % (x/A) - (y/B) = (x1/A) -(y1/B)
        xy_2 = chromosome(1:2,i+1);
        A = xy_2(1) - xy_1(1); % A = x2-x1
        B = xy_2(2) - xy_1(2); % B = y2-y1
        if A == 0
            x_curr_coe = 1;
            y_curr_coe = 0;
            c_curr = xy_1(1);
        elseif B == 0
            x_curr_coe = 0;
            y_curr_coe = 1;
            c_curr = xy_1(2);
        else
            x_curr_coe = 1/A;   % coefficient of x for path line
            y_curr_coe = -1/B;   % coefficient of y for path line
            c_curr = xy_1(1)/A - xy_1(2)/B; % c for path line
        end
        x_coe_obstacle = obstacle_xy(1,j);
        y_coe_obstacle = obstacle_xy(2,j);
        c_coe_obstacle = obstacle_c(j);
        
        if inv([x_coe_obstacle y_coe_obstacle; x_curr_coe y_curr_coe]) == Inf
            chromosome(1,end+1) = 0;
            %p=p+1;
        else
            intersection_pt = (inv([x_coe_obstacle y_coe_obstacle; x_curr_coe y_curr_coe]))*([c_coe_obstacle;c_curr]);
            
            %% Checking if intersection pt lies in obstacle
            
            if         intersection_pt(1,1)>=min(chromosome(1,i),chromosome(1,i+1))...
                    && intersection_pt(1,1)<=max(chromosome(1,i+1),chromosome(1,i))...
                    && intersection_pt(2,1)>=min(chromosome(2,i),chromosome(2,i+1))...
                    && intersection_pt(2,1)<=max(chromosome(2,i+1),chromosome(2,i))...
                    && intersection_pt(1,1)>=min(chromosome(1,i),chromosome(1,i+1))...
                    && intersection_pt(2,1)>=obstacle_c(count) && intersection_pt(2,1)<=obstacle_c(count+1)...
                    && intersection_pt(1,1)>=obstacle_c(count+2) && intersection_pt(1,1)<=obstacle_c(count+3)
                
                chromosome(1,end+1) = 1;
                
            else
                chromosome(1,end+1) = 0;
                
            end
        end
        if rem(j,4)==0
            count = count+4;
        end
    end
    %% Return value of function
    f = chromosome;
end
end

%% FITNESS FUNCTION
function fitness_chromosome = fitness(chromosome)
chromosome_local = chromosome';
Total_distance = 0;
for i = 1:1:length(chromosome_local)-1
    first_point = chromosome_local(i,:);
    second_point = chromosome_local(i+1,:);
    dist = pdist([first_point;second_point]);
    Total_distance = Total_distance + dist;
    
end
chromosome(1,end+1) = 1/Total_distance; %1/Total_distance;
fitness_chromosome =  chromosome;
end


