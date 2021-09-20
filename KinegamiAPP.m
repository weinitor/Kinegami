classdef KinegamiAPP < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        
        UIFigure                        matlab.ui.Figure
        KinegamiGEN1Label               matlab.ui.control.Label
        SelfAssignmentofJointParametersLabel  matlab.ui.control.Label
        MirroringofElbowJointsLabel     matlab.ui.control.Label
        PanelsinFinalOrigamiPrintLabel  matlab.ui.control.Label
        DXFGenerationLabel              matlab.ui.control.Label
        ElbowSplittingfortheta_mpi2Label  matlab.ui.control.Label
        OptionsSettingsLabel            matlab.ui.control.Label
        Switch                          matlab.ui.control.Switch
        Switch_2                        matlab.ui.control.Switch
        Switch_3                        matlab.ui.control.Switch
        Switch_4                        matlab.ui.control.Switch
        Switch_5                        matlab.ui.control.Switch
        
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 1150 566];
            app.UIFigure.Name = 'MATLAB App';

            % Create KinegamiGEN1Label
            app.KinegamiGEN1Label = uilabel(app.UIFigure);
            app.KinegamiGEN1Label.HorizontalAlignment = 'center';
            app.KinegamiGEN1Label.FontSize = 25;
            app.KinegamiGEN1Label.Position = [373 494 405 46];
            app.KinegamiGEN1Label.Text = 'Kinegami GEN1';

            % Create SelfAssignmentofJointParametersLabel
            app.SelfAssignmentofJointParametersLabel = uilabel(app.UIFigure);
            app.SelfAssignmentofJointParametersLabel.HorizontalAlignment = 'center';
            app.SelfAssignmentofJointParametersLabel.Position = [84 352 213 54];
            app.SelfAssignmentofJointParametersLabel.Text = 'Self-Assignment of Joint Parameters';

            % Create MirroringofElbowJointsLabel
            app.MirroringofElbowJointsLabel = uilabel(app.UIFigure);
            app.MirroringofElbowJointsLabel.HorizontalAlignment = 'center';
            app.MirroringofElbowJointsLabel.Position = [84 279 213 65];
            app.MirroringofElbowJointsLabel.Text = 'Mirroring of Elbow Joints';

            % Create PanelsinFinalOrigamiPrintLabel
            app.PanelsinFinalOrigamiPrintLabel = uilabel(app.UIFigure);
            app.PanelsinFinalOrigamiPrintLabel.HorizontalAlignment = 'center';
            app.PanelsinFinalOrigamiPrintLabel.Position = [84 213 213 49];
            app.PanelsinFinalOrigamiPrintLabel.Text = '3 Panels in Final Origami Print';

            % Create DXFGenerationLabel
            app.DXFGenerationLabel = uilabel(app.UIFigure);
            app.DXFGenerationLabel.HorizontalAlignment = 'center';
            app.DXFGenerationLabel.Position = [84 128 213 69];
            app.DXFGenerationLabel.Text = 'DXF Generation';

            % Create ElbowSplittingfortheta_mpi2Label
            app.ElbowSplittingfortheta_mpi2Label = uilabel(app.UIFigure);
            app.ElbowSplittingfortheta_mpi2Label.HorizontalAlignment = 'center';
            app.ElbowSplittingfortheta_mpi2Label.Position = [84 41 213 88];
            app.ElbowSplittingfortheta_mpi2Label.Text = 'Elbow Splitting for theta_m > pi/2';

            % Create OptionsSettingsLabel
            app.OptionsSettingsLabel = uilabel(app.UIFigure);
            app.OptionsSettingsLabel.HorizontalAlignment = 'center';
            app.OptionsSettingsLabel.FontSize = 25;
            app.OptionsSettingsLabel.Position = [84 426 346 38];
            app.OptionsSettingsLabel.Text = 'Options / Settings';

            % Create Switch
            app.Switch = uiswitch(app.UIFigure, 'slider');
            app.Switch.Position = [351 374 22 10];

            % Create Switch_2
            app.Switch_2 = uiswitch(app.UIFigure, 'slider');
            app.Switch_2.Position = [351 306 22 10];

            % Create Switch_3
            app.Switch_3 = uiswitch(app.UIFigure, 'slider');
            app.Switch_3.Position = [351 232 22 10];

            % Create Switch_4
            app.Switch_4 = uiswitch(app.UIFigure, 'slider');
            app.Switch_4.Position = [351 157 22 10];

            % Create Switch_5
            app.Switch_5 = uiswitch(app.UIFigure, 'slider');
            app.Switch_5.Position = [351 80 22 10];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = KinegamiAPP

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end