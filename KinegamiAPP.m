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
        EditFieldLabel                  matlab.ui.control.Label
        EditField                       matlab.ui.control.NumericEditField
        EditFieldLabel_2                matlab.ui.control.Label
        EditField_2                     matlab.ui.control.NumericEditField
        EditFieldLabel_3                matlab.ui.control.Label
        EditField_3                     matlab.ui.control.NumericEditField
        TextLabel                       matlab.ui.control.Label
        Button                          matlab.ui.control.Button
              
    end

    % Component initialization
    methods (Access = private)
        
%         function OpenFcn(app)
%             
%             handles = struct();
%             guidata(handles, app);
%             
%         end

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 650 566];
            app.UIFigure.Name = 'MATLAB App';

            % Create KinegamiGEN1Label
            app.KinegamiGEN1Label = uilabel(app.UIFigure);
            app.KinegamiGEN1Label.HorizontalAlignment = 'center';
            app.KinegamiGEN1Label.FontSize = 25;
            app.KinegamiGEN1Label.Position = [400 494 200 46];
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
            app.OptionsSettingsLabel.Position = [84 426 213 38];
            app.OptionsSettingsLabel.Text = 'Options / Settings';

            % Create Switch
            app.Switch = uiswitch(app.UIFigure, 'slider', ...
                'Items', {'Off', 'On'}, ...
                'ValueChangedFcn', @Switch_1Moved);
            app.Switch.Position = [351 374 22 10];

            % Create Switch_2
            app.Switch_2 = uiswitch(app.UIFigure, 'slider', ...
                'Items', {'Off', 'On'}, ...
                'ValueChangedFcn', @Switch_2Moved);
            app.Switch_2.Position = [351 306 22 10];

            % Create Switch_3
            app.Switch_3 = uiswitch(app.UIFigure, 'slider', ...
                'Items', {'Off', 'On'}, ...
                'ValueChangedFcn', @Switch_3Moved);
            app.Switch_3.Position = [351 232 22 10];

            % Create Switch_4
            app.Switch_4 = uiswitch(app.UIFigure, 'slider', ...
                'Items', {'Off', 'On'}, ...
                'ValueChangedFcn', @Switch_4Moved);
            app.Switch_4.Position = [351 157 22 10];

            % Create Switch_5
            app.Switch_5 = uiswitch(app.UIFigure, 'slider', ...
                'Items', {'Off', 'On'}, ...
                'ValueChangedFcn', @Switch_5Moved);
            app.Switch_5.Position = [351 80 22 10];
            
            % Create EditFieldLabel
            app.EditFieldLabel = uilabel(app.UIFigure);
            app.EditFieldLabel.HorizontalAlignment = 'right';
            app.EditFieldLabel.Position = [450 352 60 54];
            app.EditFieldLabel.Text = 'n  = ';

            % Create EditField
            app.EditField = uieditfield(app.UIFigure, 'numeric', ...
                'ValueChangedFcn', @(app, event) ...
                n_value_changed(app, event));
            app.EditField.Position = [515 364 40 30];
            
            % Create EditFieldLabel_2
            app.EditFieldLabel_2 = uilabel(app.UIFigure);
            app.EditFieldLabel_2.HorizontalAlignment = 'right';
            app.EditFieldLabel_2.Position = [450 279 60 65];
            app.EditFieldLabel_2.Text = '# Joints = ';
            
            % Create EditField_2
            app.EditField_2 = uieditfield(app.UIFigure, 'numeric', ...
                'ValueChangedFcn', @(app, event) ...
                nz_value_changed(app, event));
            app.EditField_2.Position = [515 296 40 30];
            
            % Create EditFieldLabel_3
            app.EditFieldLabel_3 = uilabel(app.UIFigure);
            app.EditFieldLabel_3.HorizontalAlignment = 'right';
            app.EditFieldLabel_3.Position = [450 213 60 49];
            app.EditFieldLabel_3.Text = 'r  = ';
            
            % Create EditField_3
            app.EditField_3 = uieditfield(app.UIFigure, 'numeric', ...
                'ValueChangedFcn', @(app, event) ...
                r_value_changed(app, event));
            app.EditField_3.Position = [515 222 40 30];
            
            % Create TextLabel
            app.TextLabel = uilabel(app.UIFigure);
            app.TextLabel.HorizontalAlignment = 'left';
            app.TextLabel.Position = [565 213 60 49];
            app.TextLabel.Text = '[m]';
            
            % Create Button
            app.Button = uibutton(app.UIFigure, 'push', ...
                'Callback', @Continue_Callback);
            app.Button.Position = [470 70 100 50];
            app.Button.Text = 'Continue';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
            
            % Functions which govern dynamics of GUI
            
            % Switch 1 On/Off: (Self-Assignment)
            function Switch_1Moved(Switch, event, app)
                switch Switch.Value
                    case 'Off'
                        app.Switch.Value = 'false';
                    case 'On'
                        app.Switch.selfassign = 'true';
                end
            end
            
            % Switch 2 On/Off: (Mirroring)
            function Switch_2Moved(Switch_2, event, app)
                switch Switch_2.Value
                    case 'Off'
                        app.elbow_tuck = 'false';
                    case 'On'
                        app.elbow_tuck = 'true';
                end
            end
                        
            
            % Switch 3 On/Off: (Tri-Panel)
            function Switch_3Moved(Switch_3, event, app)
                switch Switch_3.Value
                    case 'Off'
                        app.triple = 'none';
                    case 'On'
                        app.triple = 'triple';
                end
            end
            
            % Switch 4 On/Off: (DXF Generation)
            function Switch_4Moved(Switch_4, event, app)
                switch Switch_4.Value
                    case 'Off'
                        app.DXF = 'off';
                    case 'On'
                        app.DXF = 'on';
                end
            end
            
            % Switch 5 On/Off (Elbow Splitting)
            function Switch_5Moved(Switch_5, event, app)
                switch Switch_5.Value
                    case 'Off'
                        app.split = 'off';
                    case 'On'
                        app.split = 'on';
                end
            end
            
            % Text Box (n)
            function n_value_changed(txt, app)
                app.n = txt.Value;
            end            
            
            % Text Box (# joints)
            function nz_value_changed(txt, app)
                app.nz = txt.Value;
            end            
            
            % Text Box (r)
            function r_value_changed(txt, app)
                app.r = txt.Value;
            end
            
            
            % Upon Continue Button Press
            function Continue_Callback(app)
                
                disp(app.r);
                
            end

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