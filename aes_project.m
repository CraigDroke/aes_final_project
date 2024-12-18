% This script is designed to calculate annual average solar energy potential evaluator, given any
% specific location, any solar collector area, and any fixed tilt angle
% (due south)

% Authors: Craig Droke, James Galeno, Avery Mackin
% Date: November 29th, 2024

% Define a class for the interface
classdef aes_project < matlab.apps.AppBase

    % Define interface fields
    properties (Access = public)
        UIFigure             matlab.ui.Figure
        LatitudeEditField    matlab.ui.control.NumericEditField
        LongitudeEditField   matlab.ui.control.NumericEditField
        CollectorAreaField   matlab.ui.control.NumericEditField
        TiltAngleField       matlab.ui.control.NumericEditField
        SolarDataFileButton  matlab.ui.control.Button
        EvaluateButton       matlab.ui.control.Button
        ResultsTextArea      matlab.ui.control.TextArea
        SolarDataFilePath    string
    end

    % Callback functions
    methods (Access = private)

        % Evaluate the annual solar potential (calls the other function)
        function evaluatePotential(app, ~)
            % Gather user inputs from the interface
            lat = app.LatitudeEditField.Value;
            lon = app.LongitudeEditField.Value;
            area = app.CollectorAreaField.Value;
            tilt = app.TiltAngleField.Value;

            % Call the calculation function
            try
                annual_energy = calculate_solar_potential_theoretical(lat, lon, area, tilt);
                app.ResultsTextArea.Value = sprintf( ...
                    ['Results:\nLocation: Latitude %.2f, Longitude %.2f\n' ...
                    'Collector Area: %.2f m²\nTilt Angle: %.2f degrees\n' ...
                    'Annual Energy: %.2f kWh'], ...
                    lat, lon, area, tilt, annual_energy);
            % Catch the error
            catch ME
                app.ResultsTextArea.Value = ['Error: ', ME.message];
            end
        end
    end

    % App init
    methods (Access = private)
        
        % Functionto generate all of the interface components
        function createComponents(app)
            % Create the main figure
            app.UIFigure = uifigure('Position', [100, 100, 400, 300], 'Name', 'Solar Potential Evaluator');

            % Latitude input
            uilabel(app.UIFigure, 'Position', [20, 250, 80, 22], 'Text', 'Latitude:');
            app.LatitudeEditField = uieditfield(app.UIFigure, 'numeric', 'Position', [100, 250, 100, 22]);

            % Longitude input
            uilabel(app.UIFigure, 'Position', [20, 220, 80, 22], 'Text', 'Longitude:');
            app.LongitudeEditField = uieditfield(app.UIFigure, 'numeric', 'Position', [100, 220, 100, 22]);

            % Collector area input
            uilabel(app.UIFigure, 'Position', [20, 190, 120, 22], 'Text', 'Collector Area (m²):');
            app.CollectorAreaField = uieditfield(app.UIFigure, 'numeric', 'Position', [150, 190, 100, 22]);

            % Tilt angle input
            uilabel(app.UIFigure, 'Position', [20, 160, 80, 22], 'Text', 'Tilt Angle:');
            app.TiltAngleField = uieditfield(app.UIFigure, 'numeric', 'Position', [100, 160, 100, 22]);

            % Evaluate button
            app.EvaluateButton = uibutton(app.UIFigure, 'push', ...
                'Position', [20, 130, 150, 22], ...
                'Text', 'Evaluate Potential', ...
                'ButtonPushedFcn', @(~,~) app.evaluatePotential());

            % Results text area
            app.ResultsTextArea = uitextarea(app.UIFigure, ...
                'Position', [20, 20, 360, 100], ...
                'Editable', 'off', ...
                'Value', {'Results will appear here...'});
        end
    end

    % Constructor for class
    methods (Access = public)
        function app = aes_project()
            % Create and configure components
            createComponents(app);
        end
    end
end

% A function for calculating the annual solar energy at a certain latitude
% and longitude. It accepts as input the latitude in degrees, lat, the
% longitude in degrees, lon, the area in square meters, area, and the tilt
% angle, tilt, in degrees. It outputs the annual energy in kWh.
function annual_energy = calculate_solar_potential_theoretical(lat, lon, area, tilt)

    % Declaring variables
    phi_c = 0; % I think this is right, since we are always facing south
    phi_s = 0; % Sun azimuth
    A = 0;
    k = 0;
    m = 0; % Air-mass ratio
    beta = 0; % Altitude angle
    delta = 0; % Solar declination
    theta = 0; % Incidence angle
    H = 0; % Hour angle
    IB = 0; % Solar radiation
    IBC = 0; % Collected solar radiation

    day_energies = zeros(1, 365); % Vector to hold energies for each day

    % Calculating the energy (kWh) for each day in the year
    for n = 1:365

        A = 1160 + (75 * sind((360/365)*(n - 275))); % W/m^2
        k = 0.174 + (0.035 * sind((360/365)*(n - 100)));
        delta = 23.45 * sind((360/365)*(n - 81)); % Degrees

        hour_energies = zeros(1, 24);
        hours = 0:23;
        phi_values = zeros(1, 24);

        % Calculating the energy (kWh) for each hour in each day
        for hour = 0:23

            H = (12 - hour) * 15; % Hour angle in degrees
            beta = asind(cosd(lat)*cosd(delta)*cosd(H) + sind(lat)*sind(delta));
            phi_s = asind((cosd(delta)*sind(H)) / cosd(beta));
            theta = acosd(cosd(beta)*cosd(phi_s - phi_c)*sind(tilt) + sind(beta)*cosd(tilt));
            m = sqrt(((708*sind(beta))^2) + 1417) - 708*sind(beta);
            phi_values(hour + 1) = phi_s;

            if ((phi_s < 0 && hour < 12) || (phi_s > 0 && hour > 12))

                fprintf("phi_s: %.2f, day: %d, hour: %d\n", phi_s, n + 1, hour);

            end

            % OLD CODE
            % theta < 0 || theta >= 90
            %if (cosd(H) <= tand(delta) / tand(lat) || beta < 0) % Checking if phi_s is greater than 90 degrees
                    
                %hour_energies(hour + 1) = 0;
                %fprintf("phi_s > 90\n");

            %else

                IB = A * exp(-1 * k * m); % W/m^2
                IBC = IB * cosd(theta); % W/m^2

                hour_energies(hour + 1) = (IBC * area) / 1000; % kWh, assuming that IBC remains constant for each hour
                
                if (hour_energies(hour + 1) < 0) % Checking if the daily energy is negative

                    hour_energies(hour + 1) = 0;

                end

            %end

            % FOR DEBUG
            %fprintf("A: %.2f, k: %.2f, m: %.2f, phi_s: %.2f, beta: %.2f, delta: %.2f, theta: %.2f, H: %.2f, IB: %.2f, IBC: %.2f\n", A, k, m, phi_s, beta, delta, theta, H, IB, IBC);
        end

        day_energies(n) = sum(hour_energies);

        % FOR DEBUG
        %fprintf("Energy for day %d: %.2f kWh\n", n, day_energies(n));
        % Plot the data
        %if (n == 102)
            %figure;
            %bar(hours, phi_values, 'FaceColor', [0.2, 0.6, 0.8]);
            %xlabel('Hours');
            %ylabel('Phi_s');
            %grid on;
        %end

    end

    annual_energy = sum(day_energies);

    days = 1:length(day_energies);

    % Plot the data
    figure;
    bar(days, day_energies, 'FaceColor', [0.2, 0.6, 0.8]);
    xlabel('Days');
    ylabel('Daily Energy (kWh)');
    title('Annual Energy Potential Over One Year');
    grid on;

end