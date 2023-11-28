classdef bonus < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure        matlab.ui.Figure
        GridLayout      matlab.ui.container.GridLayout
        Label_3         matlab.ui.control.Label
        sourcesSpinner  matlab.ui.control.Spinner
        nsourceLabel    matlab.ui.control.Label
        Label_2         matlab.ui.control.Label
        Label           matlab.ui.control.Label
        StartButton     matlab.ui.control.Button
    end

    
    methods (Access = private)
        
        function Pmusic=music(~,dx,Index,f_c,Xr,n_source,stride)
            J=4;
            L=size(Xr,1);
            c=340;
            theta = -90:stride:90;                                       % grid
            R_x = Xr'*Xr/L;                                                % autocorrelation estimate
            a_theta = exp(-2j*pi*(Index')*(dx*sin(theta*pi/180)) / (c/f_c));                  % steer vector
            % implement eigen-decomposition and obtain the pseudo spectrum
            [EV,D] = eig(R_x);
            EVA = diag(D);
            [~,ind] = sort(EVA);
            EV = EV(:,ind);
            Un = EV(:,1:J-n_source);                                        % noise subspace (columns are eigenvectors), size: J*(J-n_source)
            Pmusic = diag(a_theta'*(Un*Un')*a_theta);
        end
        function [ss,fc]=doastft(~,X,window,overlap,nfft,fs)
            xlen=length(X);
            rown = ceil((1+nfft)/2);
            coln = 1+fix((xlen-window)/(window-overlap));
            ss = zeros(rown, coln,4);
            co=0;
            for i=1:(window-overlap):xlen-window+1
                xtw=X(i:i+window-1,:).*hamming(window);
                s=fft(xtw,nfft);
                co=co+1;
                ss(:,co,:)=s(1:rown,:);
            end
            fc=(0:rown-1)*fs/nfft;
        end
        function doa(app,X,fs,dx,n_source)
            global ax
            J = 4;                                      % number of sensors                             % number of sources
            Index = linspace(0,J-1,J);
            window = 512;
            noverlap = window/2;
            nfft = 512;
            stride=0.5;
            
            [ss,fc]=doastft(app,X,window,noverlap,nfft,fs);
            Pmusic=zeros(361,1);
            for x=(1:window/2+1)
                Xr=squeeze(ss(x,:,:));
                f_c=fc(x);
                Pmusic=Pmusic+music(app,dx,Index,f_c,Xr,n_source,stride);
            end
            theta = -90:stride:90;
            P_sm=nfft./Pmusic;
            polarplot(ax,theta*pi./180,10*log10(abs(P_sm)),'LineWidth',3);
            ax.ThetaZeroLocation = 'top';
            ax.ThetaLim=[-90,90];
            ax.ThetaTick=[-90:10:90];
            ax.FontSize = 14;
            ax.FontWeight='bold';
            ax.ThetaMinorGrid = 'on';
            P_middle = abs(P_sm(2:end-1));
            P_front = abs(P_sm(1:end-2));
            P_back = abs(P_sm(3:end));
            logic_front = (P_middle - P_front)>0;
            logic_back = (P_middle - P_back)>0;
            logic = logic_front & logic_back;
            P_middle(~logic) = min(P_middle);
            P_local = [abs(P_sm(1));P_middle;abs(P_sm(end))];
            [~,doa_Idx] = maxk(P_local,n_source);
            doa = theta(doa_Idx);
%             [~,minIdx] = min(abs(doa));
%             doa_source = doa(minIdx);
%             [~,maxIdx] = max(abs(doa));
%             interfer = doa(maxIdx);
            if n_source==1
                app.Label.Text=['The desired source DOA with MUSIC is: ',num2str(doa(1)),' deg'];
            elseif n_source==2
                app.Label.Text=['The desired source1 DOA with MUSIC is: ',num2str(doa(1)),' deg'];
                app.Label_2.Text=['The desired source2 DOA with MUSIC is: ',num2str(doa(2)),' deg'];
            elseif n_source==3
                app.Label.Text=['The desired source1 DOA with MUSIC is: ',num2str(doa(1)),' deg'];
                app.Label_2.Text=['The desired source2 DOA with MUSIC is: ',num2str(doa(2)),' deg'];
                app.Label_3.Text=['The desired source3 DOA with MUSIC is: ',num2str(doa(3)),' deg'];
            end
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            global ax start
            ax=polaraxes(app.GridLayout);
            ax.Layout.Row = [2 4];
            ax.Layout.Column= [1 7];
            ax.ThetaZeroLocation = 'top';
            ax.ThetaLim=[-90,90];
            ax.ThetaTick=[-90:10:90];
            ax.FontSize = 14;
            ax.FontWeight='bold';
            ax.ThetaMinorGrid = 'on';
            start=1;
            rec=audioDeviceReader(16000,'NumChannels',4,'Device','麦克风 (USB YDB01 Audio Effect)','SamplesPerFrame',50);
            rec();
        end

        % Button pushed function: StartButton
        function StartButtonPushed(app, event)
            global start
            if isequal(app.StartButton.Text,'Start')
                app.StartButton.Text='Stop';
                rec=audioDeviceReader(16000,'NumChannels',4,'Device','麦克风 (USB YDB01 Audio Effect)','SamplesPerFrame',8000);
                setup(rec);release(rec)
                rec.ChannelMappingSource = 'Property';
                rec.ChannelMapping = [2,1,4,3];
                while isequal(app.StartButton.Text,'Stop')&&start==1
                    fs = 16000;
                    dx=0.025;
                    data=rec();
                    if start==0
                        break
                    end
                    doa(app,data,fs,dx,app.sourcesSpinner.Value);
                    pause(0.15)
                end
            else
                app.StartButton.Text='Start';
            end
        end

        % Close request function: UIFigure
        function UIFigureCloseRequest(app, event)
            global start
            start=0;
            delete(app)
        end

        % Value changed function: sourcesSpinner
        function sourcesSpinnerValueChanged(app, event)
            value = app.sourcesSpinner.Value;
            if value==1
                app.Label_2.Text='';
                app.Label_3.Text='';
            elseif value==2
                app.Label_3.Text='';
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 901 569];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.CloseRequestFcn = createCallbackFcn(app, @UIFigureCloseRequest, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {'0.1x', '0.7x', '0.5x', 50, '1.5x', '1x', '0.7x'};
            app.GridLayout.RowHeight = {'1x', '1x', '1x', '25x'};
            app.GridLayout.ColumnSpacing = 8.79913330078125;
            app.GridLayout.Padding = [8.79913330078125 10 8.79913330078125 10];

            % Create StartButton
            app.StartButton = uibutton(app.GridLayout, 'push');
            app.StartButton.ButtonPushedFcn = createCallbackFcn(app, @StartButtonPushed, true);
            app.StartButton.FontSize = 16;
            app.StartButton.FontWeight = 'bold';
            app.StartButton.Layout.Row = [2 3];
            app.StartButton.Layout.Column = 2;
            app.StartButton.Text = 'Start';

            % Create Label
            app.Label = uilabel(app.GridLayout);
            app.Label.FontSize = 16;
            app.Label.FontWeight = 'bold';
            app.Label.Layout.Row = 1;
            app.Label.Layout.Column = [6 7];
            app.Label.Text = '';

            % Create Label_2
            app.Label_2 = uilabel(app.GridLayout);
            app.Label_2.FontSize = 16;
            app.Label_2.FontWeight = 'bold';
            app.Label_2.Layout.Row = 2;
            app.Label_2.Layout.Column = [6 7];
            app.Label_2.Text = '';

            % Create nsourceLabel
            app.nsourceLabel = uilabel(app.GridLayout);
            app.nsourceLabel.HorizontalAlignment = 'right';
            app.nsourceLabel.FontSize = 13;
            app.nsourceLabel.FontWeight = 'bold';
            app.nsourceLabel.Layout.Row = 3;
            app.nsourceLabel.Layout.Column = 3;
            app.nsourceLabel.Text = 'sources：';

            % Create sourcesSpinner
            app.sourcesSpinner = uispinner(app.GridLayout);
            app.sourcesSpinner.Limits = [1 3];
            app.sourcesSpinner.ValueChangedFcn = createCallbackFcn(app, @sourcesSpinnerValueChanged, true);
            app.sourcesSpinner.FontSize = 13;
            app.sourcesSpinner.FontWeight = 'bold';
            app.sourcesSpinner.Layout.Row = 3;
            app.sourcesSpinner.Layout.Column = 4;
            app.sourcesSpinner.Value = 2;

            % Create Label_3
            app.Label_3 = uilabel(app.GridLayout);
            app.Label_3.FontSize = 16;
            app.Label_3.FontWeight = 'bold';
            app.Label_3.Layout.Row = 3;
            app.Label_3.Layout.Column = [6 7];
            app.Label_3.Text = '';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = bonus

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

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