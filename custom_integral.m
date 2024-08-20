
function integral_result = custom_integral(func_handle, lower_limit, upper_limit, num_intervals)
    % 使用梯形法则计算积分

    % 计算步长
    step_size = (upper_limit - lower_limit) / num_intervals;

    % 初始化积分结果
    integral_result = 0;

    % 计算积分
    for i = 1:num_intervals
        % 计算当前区间的左右端点
        x0 = lower_limit + (i - 1) * step_size;
        x1 = lower_limit + i * step_size;

        % 计算当前区间的梯形面积并加到积分结果中
        integral_result = integral_result + (func_handle(x0) + func_handle(x1)) / 2 * step_size;
    end
end
