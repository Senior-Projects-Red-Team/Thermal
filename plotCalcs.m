function[time, mean_interior_temp, mean_inner_wall, mean_outer_wall, total_power, derivative, R_conv, R_bridge, R_foam] = plotCalcs(data_table, table_title, correction, plotBool)

    time = linspace(0,height(data_table)-1, height(data_table)).*10;
    
    mean_interior_temp = zeros(height(data_table),1);
    mean_inner_wall = mean_interior_temp;
    mean_outer_wall = mean_interior_temp;
    total_power = mean_interior_temp;
    
    for i = 1:height(data_table)
        mean_interior_temp(i) = (data_table{i,"TMP1"}+data_table{i,"TMP2"})/2;
        mean_inner_wall(i) = (data_table{i,"TC2"}+data_table{i,"TC0"})/2;
        mean_outer_wall(i) = (data_table{i,"TC1"}+data_table{i,"TC3"})/2;
        total_power(i) = 16*(data_table{i,"Pin10"} + 2*data_table{i,"Pin11"} + data_table{i,"Pin6"} + data_table{i,"Pin9"})/100;
    end
    
    mean_inner_wall = ((mean_inner_wall + 273.15).*correction)-273.15;
    mean_outer_wall = ((mean_outer_wall + 273.15).*correction)-273.15;
    
    % Removing jump discontinuities.
    % for i = 10:length(mean_interior_temp)-5
    %     mean_delta = mean(mean_inner_wall(i+1:i+3))-mean(mean_inner_wall(i-2:i));
    %     if abs(mean_delta) > 4
    %         mean_inner_wall(i+1:end) = mean_inner_wall(i+1:end) - mean_delta;
    %     end
    %     mean_delta = mean(mean_outer_wall(i+1:i+3))-mean(mean_outer_wall(i-2:i));
    %     if abs(mean_delta) > 4
    %         mean_outer_wall(i+1:end) = mean_outer_wall(i+1:end) - mean_delta;
    %     end
    % end

    if plotBool
        figure()
        hold on
        title("Mean temperatures over time " + table_title)
        plot(time,mean_interior_temp)
        plot(time,mean_inner_wall)
        plot(time,mean_outer_wall)
        legend(["Interior","Inner Wall","Outer Wall"], Location="southeast")
        ylim([-50,35])
        ylabel("Temp ^oC")
        xlabel("Time (sec)")
    
    
        figure()
        hold on
        title("Total Power usage over time " + table_title)
        plot(time, smoothdata(total_power, 'movmean', 40))
        xlabel("Time (sec)")
    end
    
    derivative = (total_power(2:end) - total_power(1:end-1))./10;
    
    if plotBool
        figure()
        hold on
        title("Change in Total Power use over time (W/s) " + table_title)
        plot(time(2:end), smoothdata(derivative, "movmean", 50))
        yline(0)
        xlabel("Time (sec)")
    end
    
    R_conv = abs(mean_interior_temp - mean_inner_wall)./total_power;
    
    if plotBool
        figure()
        hold on
        title("Estimate for R_{conv} over time " + table_title)
        plot(time, R_conv)
        yline(0.497)
        legend(["From Data", "Estimate"])
        xlabel("Time (sec)")
    end
    
    R_bridge = abs(mean_outer_wall - mean_inner_wall)./total_power;
    
    if plotBool
        figure()
        hold on
        title("Estimate for R_{bridge} over time " + table_title)
        plot(time, R_bridge)
        yline(0.66)
        legend(["From Data", "Estimate"])
        xlabel("Time (sec)")
    end
    
    R_foam = abs(mean_outer_wall - -78.5)./total_power;
    
    if plotBool
    figure()
    hold on
    title("Estimate for R_{foam} over time " + table_title)
    plot(time, R_foam)
    yline(1.8)
    legend(["From Data", "Estimate"])
    end

end