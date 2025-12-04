from crewai import Agent, Crew, Process, Task
from crewai.project import CrewBase, agent, crew, task
from Latest_Rt.tools import (
    extract_literature_rt,
    predict_retention_time,
    compare_rt_values,
    batch_extract_literature_rt,
    batch_predict_retention_times,
    batch_compare_rt_values,
    generate_summary_report,
    create_performance_visualization,
    export_results_to_csv,
    generate_compound_class_report
)


@CrewBase
class LatestRt():
    """LatestRt crew"""

    agents_config = 'config/agents.yaml'
    tasks_config = 'config/tasks.yaml'

    @agent
    def literature_agent(self) -> Agent:
        return Agent(
            # type: ignore[index]
            config=self.agents_config['literature_agent'],
            tools=[extract_literature_rt, batch_extract_literature_rt],
            verbose=True
        )

    @agent
    def prediction_agent(self) -> Agent:
        return Agent(
            # type: ignore[index]
            config=self.agents_config['prediction_agent'],
            tools=[predict_retention_time, batch_predict_retention_times],
            verbose=True
        )

    @agent
    def validation_agent(self) -> Agent:
        return Agent(
            # type: ignore[index]
            config=self.agents_config['validation_agent'],
            tools=[compare_rt_values, batch_compare_rt_values, generate_summary_report,
                   create_performance_visualization, export_results_to_csv,
                   generate_compound_class_report],
            verbose=True
        )

    @agent
    def coordinator_agent(self) -> Agent:
        return Agent(
            # type: ignore[index]
            config=self.agents_config['coordinator_agent'],
            verbose=True
        )

    @task
    def extract_literature_rt_task(self) -> Task:
        return Task(
            # type: ignore[index]
            config=self.tasks_config['extract_literature_rt'],
        )

    @task
    def predict_retention_time_task(self) -> Task:
        return Task(
            # type: ignore[index]
            config=self.tasks_config['predict_retention_time'],
            output_file='report.md'
        )

    @task
    def compare_rt_values_task(self) -> Task:
        return Task(
            # type: ignore[index]
            config=self.tasks_config['compare_rt_values'],
            output_file='report.md'
        )

    @crew
    def crew(self) -> Crew:

        return Crew(
            agents=self.agents,  # Automatically created by the @agent decorator
            tasks=self.tasks,  # Automatically created by the @task decorator
            process=Process.sequential,
            verbose=True,
            max_rpm=10,  # Limit to 10 requests per minute to avoid rate limits
            # process=Process.hierarchical, # In case you wanna use that instead https://docs.crewai.com/how-to/Hierarchical/
        )
